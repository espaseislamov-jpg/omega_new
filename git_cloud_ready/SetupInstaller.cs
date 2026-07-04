using System;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Reflection;
using System.Windows.Forms;

internal static class SetupInstaller
{
    private const string AppName = "Omega v1";
    private const string AppFolderName = "Omega_v1";
    private const string PayloadResourceName = "Omega_v1_payload.zip";

    [STAThread]
    private static void Main()
    {
        Application.EnableVisualStyles();
        Application.SetCompatibleTextRenderingDefault(false);

        string installDir = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
            "Programs",
            AppFolderName
        );

        DialogResult confirm = MessageBox.Show(
            "Install " + AppName + " to:\n\n" + installDir + "\n\nContinue?",
            AppName + " Setup",
            MessageBoxButtons.OKCancel,
            MessageBoxIcon.Information
        );
        if (confirm != DialogResult.OK)
        {
            return;
        }

        try
        {
            Install(installDir);

            DialogResult launch = MessageBox.Show(
                AppName + " installed successfully.\n\nLaunch now?",
                AppName + " Setup",
                MessageBoxButtons.YesNo,
                MessageBoxIcon.Information
            );

            if (launch == DialogResult.Yes)
            {
                string exePath = Path.Combine(installDir, "Omega_v1.exe");
                if (File.Exists(exePath))
                {
                    Process.Start(exePath);
                }
            }
        }
        catch (Exception ex)
        {
            MessageBox.Show(
                ex.Message,
                AppName + " Setup",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error
            );
        }
    }

    private static void Install(string installDir)
    {
        string tempDir = Path.Combine(Path.GetTempPath(), AppFolderName + "_setup");
        if (Directory.Exists(tempDir))
        {
            Directory.Delete(tempDir, true);
        }
        Directory.CreateDirectory(tempDir);

        if (Directory.Exists(installDir))
        {
            Directory.Delete(installDir, true);
        }
        Directory.CreateDirectory(installDir);

        using (Stream payloadStream = Assembly.GetExecutingAssembly().GetManifestResourceStream(PayloadResourceName))
        {
            if (payloadStream == null)
            {
                throw new InvalidOperationException("Embedded payload not found.");
            }

            using (var archive = new ZipArchive(payloadStream, ZipArchiveMode.Read))
            {
                foreach (ZipArchiveEntry entry in archive.Entries)
                {
                    string destinationPath = Path.Combine(tempDir, entry.FullName);
                    string destinationDir = Path.GetDirectoryName(destinationPath);
                    if (!string.IsNullOrEmpty(destinationDir))
                    {
                        Directory.CreateDirectory(destinationDir);
                    }

                    if (string.IsNullOrEmpty(entry.Name))
                    {
                        continue;
                    }

                    entry.ExtractToFile(destinationPath, true);
                }
            }
        }

        CopyDirectory(tempDir, installDir);
        CreateShortcuts(installDir);
        Directory.Delete(tempDir, true);
    }

    private static void CopyDirectory(string sourceDir, string destinationDir)
    {
        foreach (string directory in Directory.GetDirectories(sourceDir, "*", SearchOption.AllDirectories))
        {
            string targetDir = directory.Replace(sourceDir, destinationDir);
            Directory.CreateDirectory(targetDir);
        }

        foreach (string file in Directory.GetFiles(sourceDir, "*", SearchOption.AllDirectories))
        {
            string targetFile = file.Replace(sourceDir, destinationDir);
            string targetDir = Path.GetDirectoryName(targetFile);
            if (!string.IsNullOrEmpty(targetDir))
            {
                Directory.CreateDirectory(targetDir);
            }
            File.Copy(file, targetFile, true);
        }
    }

    private static void CreateShortcuts(string installDir)
    {
        string exePath = Path.Combine(installDir, "Omega_v1.exe");
        string uninstallScriptPath = Path.Combine(installDir, "uninstall_release.ps1");
        string desktopShortcut = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.DesktopDirectory),
            AppName + ".lnk"
        );
        string startMenuDir = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.StartMenu),
            "Programs",
            AppName
        );
        Directory.CreateDirectory(startMenuDir);

        string startMenuShortcut = Path.Combine(startMenuDir, AppName + ".lnk");
        string uninstallShortcut = Path.Combine(startMenuDir, "Uninstall " + AppName + ".lnk");
        string powershellPath = Path.Combine(
            Environment.GetFolderPath(Environment.SpecialFolder.Windows),
            "System32",
            "WindowsPowerShell",
            "v1.0",
            "powershell.exe"
        );

        CreateShortcut(desktopShortcut, exePath, installDir, exePath + ",0", null);
        CreateShortcut(startMenuShortcut, exePath, installDir, exePath + ",0", null);
        CreateShortcut(
            uninstallShortcut,
            powershellPath,
            installDir,
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.System), "shell32.dll") + ",131",
            "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File \"" + uninstallScriptPath + "\""
        );
    }

    private static void CreateShortcut(string shortcutPath, string targetPath, string workingDirectory, string iconLocation, string arguments)
    {
        Type shellType = Type.GetTypeFromProgID("WScript.Shell");
        object shell = Activator.CreateInstance(shellType);
        object shortcut = shellType.InvokeMember("CreateShortcut", BindingFlags.InvokeMethod, null, shell, new object[] { shortcutPath });
        Type shortcutType = shortcut.GetType();
        shortcutType.InvokeMember("TargetPath", BindingFlags.SetProperty, null, shortcut, new object[] { targetPath });
        shortcutType.InvokeMember("WorkingDirectory", BindingFlags.SetProperty, null, shortcut, new object[] { workingDirectory });
        shortcutType.InvokeMember("IconLocation", BindingFlags.SetProperty, null, shortcut, new object[] { iconLocation });
        if (!string.IsNullOrWhiteSpace(arguments))
        {
            shortcutType.InvokeMember("Arguments", BindingFlags.SetProperty, null, shortcut, new object[] { arguments });
        }
        shortcutType.InvokeMember("Save", BindingFlags.InvokeMethod, null, shortcut, null);
    }
}
