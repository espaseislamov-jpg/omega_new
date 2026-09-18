#define MyAppName "Omega V3.0"
#define MyAppVersion "3.0"
#define MyAppPublisher "Omega"
#define MyAppExeName "Omega.exe"
#define MyRuntimeDir "_runtime_v3_0_1_open"

[Setup]
AppId={{8F342910-9F62-46E8-88D6-0A3E2C0B2002}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={localappdata}\Programs\Omega V3.0
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\installer_output
OutputBaseFilename=Omega_V3.0_Setup
VersionInfoVersion=3.0.1.0
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
UninstallDisplayIcon={app}\{#MyAppExeName}
; CI verification extracts this exact installer into an isolated directory
; without replacing the developer's existing uninstall registration or shortcuts.
Uninstallable=not IsVerificationInstall

[Languages]
Name: "russian"; MessagesFile: "compiler:Languages\Russian.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Создать ярлык на рабочем столе"; GroupDescription: "Дополнительные значки:"; Flags: unchecked

[Files]
; Editable profiles, training data and manual work are created by the app,
; never shipped from a developer's working folder or overwritten on update.
Source: "..\dist\Omega_V3.0\Omega.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\dist\Omega_V3.0\{#MyRuntimeDir}\*"; DestDir: "{app}\{#MyRuntimeDir}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Check: not IsVerificationInstall
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon; Check: not IsVerificationInstall

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Запустить {#MyAppName}"; Flags: nowait postinstall skipifsilent

[Code]
function IsVerificationInstall: Boolean;
begin
  Result := ExpandConstant('{param:VERIFYINSTALL|0}') = '1';
end;
