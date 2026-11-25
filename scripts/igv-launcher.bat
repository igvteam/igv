setlocal

::%~dp0 is the directory containing this file.  It is quoted rather than converted to a Windows 8.3 "short name"
::as igv.bat does, as the short name may not exist for all path components, particularly on drives other than C:,
::which lets spaces into the path (e.g. "Program Files").

if exist "%~dp0jdk-21" (
  echo "Using bundled JDK."
  set "JAVA_HOME=%~dp0jdk-21"
  set JAVA_CMD="%~dp0jdk-21\bin\javaw"
) else (
  echo "Using system JDK."
  set JAVA_CMD=java
)

::-Xmx8g indicates 8 gb of memory.
::To adjust this (or other Java options), edit the "%USERPROFILE%\.igv\java_arguments"
::file.  For more info, see the README at
::https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt
if exist "%USERPROFILE%\.igv\java_arguments" (
  start "" %JAVA_CMD% -showversion -Xmx8g @"%~dp0igv.args" @"%USERPROFILE%\.igv\java_arguments" -cp "%~dp0lib\*" org.igv.ui.Main %*
) else (
  start "" %JAVA_CMD% -showversion -Xmx8g @"%~dp0igv.args" -cp "%~dp0lib\*" org.igv.ui.Main %*
)
