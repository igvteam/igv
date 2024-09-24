setlocal

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
  start %JAVA_CMD% -showversion --module-path=lib -Xmx8g @igv.args @"%USERPROFILE%\.igv\java_arguments" --module=org.igv/org.broad.igv.ui.Main
) else (
  start %JAVA_CMD% -showversion --module-path=lib -Xmx8g @igv.args --module=org.igv/org.broad.igv.ui.Main
)
