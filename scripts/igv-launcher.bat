setlocal

::This launcher is started by the shortcuts the installer creates, in the installation directory.  It uses
::relative paths, as the Windows 8.3 "short name" igv.bat relies on may not exist for all path components,
::particularly on drives other than C:, which lets spaces into the path (e.g. "Program Files").

if exist jdk-21 (
  echo "Using bundled JDK."
  set JAVA_HOME=jdk-21
  set JAVA_CMD=jdk-21\bin\javaw
) else (
  echo "Using system JDK."
  set JAVA_CMD=java
)

::-Xmx8g indicates 8 gb of memory.
::To adjust this (or other Java options), edit the "%USERPROFILE%\.igv\java_arguments"
::file.  For more info, see the README at
::https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt
if exist "%USERPROFILE%\.igv\java_arguments" (
  start %JAVA_CMD% -showversion -Xmx8g @igv.args @"%USERPROFILE%\.igv\java_arguments" -cp "lib\*" org.igv.ui.Main %*
) else (
  start %JAVA_CMD% -showversion -Xmx8g @igv.args -cp "lib\*" org.igv.ui.Main %*
)
