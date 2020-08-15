setlocal
::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

if exist %BatchPath%\jdk-11 (
  echo "Using bundled JDK."
  set JAVA_HOME=%BatchPath%\jdk-11
  set JAVA_CMD=%BatchPath%\jdk-11\bin\javaw
) else (
  echo "Using system JDK."
  set JAVA_CMD=java
)

::-Xmx4g indicates 4 gb of memory.
::To adjust this (or other Java options), edit the "%USERPROFILE%\.igv\java_arguments" 
::file.  For more info, see the README at 
::https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt 
if exist "%USERPROFILE%\.igv\java_arguments" (
  start %JAVA_CMD% -showversion --module-path=lib -Xmx4g -Dproduction=true @igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true @"%USERPROFILE%\.igv\java_arguments" --module=org.igv/org.broad.igv.ui.Main 
) else (
  start %JAVA_CMD% -showversion --module-path=lib -Xmx4g -Dproduction=true @igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true --module=org.igv/org.broad.igv.ui.Main 
)
