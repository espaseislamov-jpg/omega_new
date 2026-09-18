"""Immutable signal snapshots and operator reviews with optional automatic grades."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import uuid
import numpy as np
from omega_core.io import get_runtime_app_dir
from omega_core import runtime

SCHEMA=1
POINTS=96
HALF_WINDOW=.08
ISSUES={
    'left':'Левая граница', 'right':'Правая граница',
    'shoulder':'Захвачено плечо / соседний пик', 'missing':'Пропущен пик',
    'identity':'Назначена другая кислота', 'separation':'Пики не разделены',
    'other':'Другая проблема',
}


def workspace():return get_runtime_app_dir()/'judge_workspace'


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    temporary.write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    temporary.replace(path)


def fingerprint(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def finite(value,default):
    try:value=float(value)
    except (ValueError,TypeError):return default
    return value if np.isfinite(value) else default


def features(signal,rows):
    """Fixed physical-time windows of raw signal, corrected signal and marks.

    No old judge score, Omega result, filename, manual answer or human grade.
    """
    x,raw,baseline=signal.T
    if len(x)<2 or not np.isfinite(signal).all() or np.any(np.diff(x)<=0):
        raise ValueError('Некорректный исходный сигнал для нейросудьи.')
    output=[]
    for row in rows:
        found=finite(row.get('found_rt'),np.nan)
        expected=finite(row.get('corrected_target_rt'),finite(row.get('expected_rt'),np.nan))
        center=found if np.isfinite(found) else expected
        if not np.isfinite(center):center=float(x[0])
        grid=np.linspace(center-HALF_WINDOW,center+HALF_WINDOW,POINTS)
        r=np.interp(grid,x,raw);b=np.interp(grid,x,baseline)
        corrected=r-b
        scale=max(float(np.max(abs(corrected))),float(np.ptp(r)),1e-8)
        a=finite(row.get('integration_start_x'),np.nan);z=finite(row.get('integration_end_x'),np.nan)
        mask=((grid>=a)&(grid<=z)).astype(float)
        apex=np.exp(-.5*((grid-found)/.0015)**2) if np.isfinite(found) else np.zeros(POINTS)
        channels=np.array([(r-np.median(b))/scale,corrected/scale,mask,apex])
        context=[float(np.isfinite(found)),float(np.isfinite(a) and np.isfinite(z) and a<z),
                 finite((a-center)/HALF_WINDOW,0),finite((z-center)/HALF_WINDOW,0),
                 finite((found-expected)/HALF_WINDOW,0),float((grid>=x[0]).all() and (grid<=x[-1]).all())]
        output.append(np.r_[channels.ravel(),np.clip(context,-4,4)])
    return np.asarray(output,dtype=np.float32)


def capture(batch,source,profile):
    processed=batch['processed_df']
    signal=processed[['x_corrected','y','baseline']].to_numpy(dtype=np.float64,copy=True)
    frame=batch['matched_targets_df'].sort_values('order_index')
    rows=clean(frame.to_dict('records'))
    source=Path(source)
    date=re.search(r'(?<!\d)(\d{8})(?!\d)',source.stem)
    source_hash=fingerprint(source)
    sample=str(batch['sample_name'])
    metadata=clean(dict(schema=SCHEMA,source_name=source.name,source_sha256=source_hash,
        sample=sample,sample_identity=re.sub(r'^O\d+_','',sample).removesuffix('.D').casefold(),
        batch_group=date.group(1) if date else source_hash,
        profile=profile,engine_sha256=runtime.code_digest(),baseline_mode=batch.get('baseline_mode'),
        omega=batch.get('omega_report'),rows=rows,codes=frame.code.astype(str).tolist()))
    digest=hashlib.sha256(json.dumps(metadata,sort_keys=True,ensure_ascii=False).encode())
    digest.update(signal.astype('<f8').tobytes())
    metadata['state_id']=digest.hexdigest()
    return dict(metadata=metadata,signal=signal,features=features(signal,rows))


def default_labels(codes):
    return {code:dict(verdict='good',issues=[],note='',serious=False) for code in codes}


def automatic_rating(peak_labels):
    """Each issue counts once; a serious peak contributes at least four points."""
    points=0
    for item in peak_labels.values():
        if item.get('verdict')!='bad':continue
        count=max(1,len(set(item.get('issues',[]))))
        points+=max(4,count) if item.get('serious',False) else count
    return max(0,5-(points+1)//2),points


def validate_review(snapshot,rating,peak_labels,note=''):
    if rating is not None and (type(rating) is not int or not 0<=rating<=5):
        raise ValueError('Выберите оценку от 0 до 5 или автоматическую оценку.')
    known=set(snapshot['metadata']['codes'])
    result={}
    for code,item in peak_labels.items():
        if code not in known:raise ValueError('В разметке найден неизвестный пик.')
        verdict=item.get('verdict','unreviewed');issues=list(dict.fromkeys(item.get('issues',[])))
        if verdict not in {'unreviewed','good','bad','unsure'} or set(issues)-set(ISSUES):
            raise ValueError('Некорректная отметка пика.')
        if issues and verdict!='bad':raise ValueError('Пик с замечаниями должен иметь отметку «Есть замечания».')
        if verdict=='bad' and not issues:raise ValueError('Для проблемного пика отметьте причину, хотя бы «Другая проблема».')
        serious=item.get('serious',False)
        if type(serious) is not bool or (serious and verdict!='bad'):
            raise ValueError('Серьёзность можно отметить только у проблемного пика.')
        result[code]=dict(verdict=verdict,issues=issues,note=str(item.get('note','')).strip(),serious=serious)
    auto,points=automatic_rating(result)
    return dict(rating=auto if rating is None else rating,manual_rating=rating,
                rating_source='automatic' if rating is None else 'manual',
                automatic_rating=auto,remark_points=points,rating_rule='issue_count_serious4_v1',
                peaks=result,note=str(note).strip())


def save_review(snapshot,rating,peak_labels,note='',root=None,prediction_seen=None):
    root=Path(root) if root is not None else workspace()
    label=validate_review(snapshot,rating,peak_labels,note)
    folder=root/'examples'/snapshot['metadata']['state_id'];folder.mkdir(parents=True,exist_ok=True)
    if not (folder/'signal.npz').exists():
        temporary=folder/(uuid.uuid4().hex+'.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream,signal=snapshot['signal'],features=snapshot['features'])
        temporary.replace(folder/'signal.npz')
    atomic_json(folder/'snapshot.json',snapshot['metadata'])
    label.update(schema=SCHEMA,state_id=snapshot['metadata']['state_id'],
        review_id=uuid.uuid4().hex,created_at=datetime.now(timezone.utc).isoformat(),prediction_seen=prediction_seen)
    atomic_json(folder/'history'/(label['review_id']+'.json'),label)
    atomic_json(folder/'review.json',label)
    return folder


def existing_review(snapshot,root=None):
    root=Path(root) if root is not None else workspace()
    path=root/'examples'/snapshot['metadata']['state_id']/'review.json'
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else None


def load_examples(root=None):
    root=Path(root) if root is not None else workspace()
    result=[]
    for path in sorted((root/'examples').glob('*/review.json')):
        meta=json.loads((path.parent/'snapshot.json').read_text(encoding='utf-8'))
        review=json.loads(path.read_text(encoding='utf-8'))
        if meta['schema']!=SCHEMA or review['state_id']!=meta['state_id']:
            raise ValueError('Версия или идентификатор сохранённой разметки не совпадает.')
        with np.load(path.parent/'signal.npz',allow_pickle=False) as arrays:
            x=arrays['features'].copy()
        validate_review({'metadata':meta},review['rating'],review['peaks'],review.get('note',''))
        if not np.isfinite(x).all():raise ValueError('Повреждены признаки сохранённого примера.')
        result.append(dict(metadata=meta,review=review,x=x))
    from .quiz_data import load_imported_examples
    return result + load_imported_examples(root)


def split_groups(examples):
    """Union batches sharing sample identities or a source file; never leak variants."""
    parents=list(range(len(examples)))
    def find(i):
        while parents[i]!=i:
            parents[i]=parents[parents[i]];i=parents[i]
        return i
    owners={}
    for i,e in enumerate(examples):
        m=e['metadata']
        for key in [('batch',m['batch_group']),('sample',m['sample_identity']),('source',m['source_sha256'])]:
            if key in owners:parents[find(i)]=find(owners[key])
            else:owners[key]=i
    return np.array([find(i) for i in range(len(examples))])
