"""Portable, explicit partial peak reviews; no automatic training or transmission."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import uuid
import zipfile
import numpy as np
from . import data


def outbox(root=None):
    return (Path(root) if root is not None else data.workspace())/'quiz_outbox'


def choose_questions(app, count=6, root=None):
    root = Path(root) if root is not None else data.workspace()
    pool = []
    for path in (root/'examples').glob('*/snapshot.json'):
        metadata = json.loads(path.read_text(encoding='utf-8'))
        for row in metadata['rows']:
            pool.append((path, metadata, row['code']))
    rng = random.SystemRandom()
    rng.shuffle(pool)
    questions, used_codes, loaded = [], set(), {}
    for path, metadata, code in pool:
        if code in used_codes:
            continue
        if path not in loaded:
            with np.load(path.parent/'signal.npz', allow_pickle=False) as arrays:
                signal = arrays['signal'].copy()
            loaded[path] = dict(metadata=metadata, signal=signal, features=data.features(signal,metadata['rows']))
        questions.append(dict(snapshot=loaded[path], code=code))
        used_codes.add(code)
        if len(questions)==count:
            return questions
    # Colleagues can contribute from the current batch without an installed model.
    batches = [b for b in app.loaded_batches if b.get('processed_df') is not None]
    rng.shuffle(batches)
    for batch in batches:
        snapshot = data.capture(batch, app.current_file, app.active_instrument_profile)
        codes = list(snapshot['metadata']['codes'])
        rng.shuffle(codes)
        for code in codes:
            if code not in used_codes:
                questions.append(dict(snapshot=snapshot, code=code))
                used_codes.add(code)
                if len(questions)==count:
                    return questions
    return questions


def validate_snapshot(metadata, signal):
    if not isinstance(metadata, dict) or metadata.get('schema') != data.SCHEMA:
        raise ValueError('Неизвестный формат снимка.')
    signal = np.asarray(signal, dtype=np.float64)
    if signal.ndim != 2 or signal.shape[1] != 3 or not 2 <= len(signal) <= 250000:
        raise ValueError('Некорректный размер сигнала.')
    codes = metadata['codes']
    if not 1 <= len(codes) <= 100 or len(set(codes))!=len(codes) or codes != [r['code'] for r in metadata['rows']]:
        raise ValueError('Не совпадают кислоты снимка.')
    original = {k:v for k,v in metadata.items() if k!='state_id'}
    digest = hashlib.sha256(json.dumps(original,sort_keys=True,ensure_ascii=False).encode())
    digest.update(signal.astype('<f8').tobytes())
    if metadata['state_id'] != digest.hexdigest():
        raise ValueError('Контрольная сумма снимка не совпала.')
    return dict(metadata=metadata, signal=signal, features=data.features(signal,metadata['rows']))


def validate_package(payload):
    if payload.get('format') != 'omega-peak-quiz' or payload.get('version') != 1:
        raise ValueError('Неизвестный формат квиза.')
    if not isinstance(payload.get('session_id'),str) or len(payload['session_id'])!=32 or any(c not in '0123456789abcdef' for c in payload['session_id']):
        raise ValueError('Некорректный идентификатор квиза.')
    if not 1 <= len(payload['answers']) <= 6 or not 1 <= len(payload['snapshots']) <= 6:
        raise ValueError('В квизе должно быть от 1 до 6 ответов.')
    snapshots = {}
    for entry in payload['snapshots']:
        snapshot = validate_snapshot(entry['metadata'],entry['signal'])
        state = snapshot['metadata']['state_id']
        if state in snapshots:
            raise ValueError('Повторный снимок.')
        snapshots[state] = snapshot
    seen = set()
    for answer in payload['answers']:
        key = (answer['state_id'],answer['code'])
        if key in seen or key[0] not in snapshots:
            raise ValueError('Повторный ответ или неизвестный снимок.')
        seen.add(key)
        data.validate_review(snapshots[key[0]], None, {key[1]:answer['label']})
        grade = answer.get('peak_rating')
        if grade is not None and (type(grade) is not int or not 0<=grade<=5):
            raise ValueError('Оценка пика должна быть от 0 до 5.')
    return snapshots


def save_answers(answers, operator='', root=None):
    if not answers:
        return None
    snapshots, records = {}, []
    for item in answers:
        snapshot = item['snapshot']
        state = snapshot['metadata']['state_id']
        snapshots[state] = dict(metadata=snapshot['metadata'], signal=snapshot['signal'].tolist())
        records.append(dict(state_id=state, code=item['code'], label=item['label'], peak_rating=item.get('peak_rating')))
    payload = dict(format='omega-peak-quiz', version=1, session_id=uuid.uuid4().hex,
                   created_at=datetime.now(timezone.utc).isoformat(), operator=str(operator),
                   snapshots=list(snapshots.values()), answers=records, prediction_seen=False)
    validate_package(payload)
    directory = outbox(root)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory/f"quiz_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{payload['session_id'][:8]}.omega-quiz"
    temporary = path.with_suffix('.tmp')
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('quiz.json',json.dumps(data.clean(payload),ensure_ascii=False,allow_nan=False))
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def import_package(path, root=None):
    root = Path(root) if root is not None else data.workspace()
    path = Path(path)
    if path.stat().st_size > 32_000_000:
        raise ValueError('Слишком большой файл квиза.')
    with zipfile.ZipFile(path) as archive:
        if archive.namelist()!=['quiz.json'] or archive.getinfo('quiz.json').file_size>100_000_000:
            raise ValueError('Некорректное содержимое файла квиза.')
        payload = json.loads(archive.read('quiz.json'))
    validate_package(payload)
    destination = root/'quiz_imports'/(payload['session_id']+'.json')
    if destination.exists():
        if json.loads(destination.read_text(encoding='utf-8')) != payload:
            raise ValueError('Конфликт: такой номер квиза уже импортирован с другим содержимым.')
        return 0
    data.atomic_json(destination, payload)
    return len(payload['answers'])


def load_imported_examples(root):
    examples = []
    for path in (Path(root)/'quiz_imports').glob('*.json'):
        payload = json.loads(path.read_text(encoding='utf-8'))
        snapshots = validate_package(payload)
        for state, snapshot in snapshots.items():
            answers = [a for a in payload['answers'] if a['state_id']==state]
            labels = {a['code']:a['label'] for a in answers}
            # A few reviewed peaks do not constitute a whole-sample grade.
            review = dict(rating=None, rating_source='partial_quiz', peaks=labels,
                          review_id=payload['session_id'], prediction_seen=None,
                          peak_ratings={a['code']:a.get('peak_rating') for a in answers})
            examples.append(dict(metadata=snapshot['metadata'], review=review, x=snapshot['features']))
    return examples
