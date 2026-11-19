# Workaround for CVE-2026-114938 (Java4Shell)

[CVE-2026-114938](https://www.cve.org/CVERecord?id=CVE-2026-114938) is a vulnerability in Java 17 and later, 
including the Java runtime bundled with current and recent IGV releases.  It will be fixed in an upcoming IGV 
release.  Until you can upgrade, the workaround below turns off the affected Java feature.  It does not change 
how IGV works.

The vulnerability only affects computers with certain Intel and AMD processors.  Macs with Apple Silicon 
(M1 and later) are not affected and do not need the workaround.

## If you already have a java_arguments file

Do not replace it.  Open it in a text editor, add this as a new line, and restart IGV:

```
-XX:-UseBASE64Intrinsics
```

## If you do not have a java_arguments file

Save the file [java_arguments](https://raw.githubusercontent.com/igvteam/igv/main/notes/java4shell/java_arguments)
in a folder named `.igv` in your home folder, then restart IGV.

| | Location of the file |
|---|---|
| Windows 10 and 11 | `C:\Users\<your user name>\.igv\java_arguments` |
| macOS | `/Users/<your user name>/.igv/java_arguments` |
| Linux | `/home/<your user name>/.igv/java_arguments` |

Notes

* The file name must be exactly `java_arguments`, with no extension.  It will not work if it is saved as 
  `java_arguments.txt`.
* The `.igv` folder is **not** the `igv` folder that holds your preferences and genomes.  It is a separate 
  folder, the name starts with a period, and it usually does not exist until you create it.

The commands below create the folder and download the file to the correct place with the correct name.

### Windows

Open **Command Prompt** (press the Windows key, type `cmd`, press Enter) and paste in these two lines, pressing 
Enter after each.  Use Command Prompt, not PowerShell.

```
if not exist "%USERPROFILE%\.igv" mkdir "%USERPROFILE%\.igv"
curl -L -o "%USERPROFILE%\.igv\java_arguments" https://raw.githubusercontent.com/igvteam/igv/main/notes/java4shell/java_arguments
```

### macOS and Linux

Open a terminal (on macOS, the **Terminal** application) and paste in these two lines, pressing Enter after each.

```
mkdir -p "$HOME/.igv"
curl -L -o "$HOME/.igv/java_arguments" https://raw.githubusercontent.com/igvteam/igv/main/notes/java4shell/java_arguments
```

## Checking the result

The file should contain the line `-XX:-UseBASE64Intrinsics`.  To check, run this in the same window:

* Windows: `type "%USERPROFILE%\.igv\java_arguments"`
* macOS and Linux: `cat "$HOME/.igv/java_arguments"`
