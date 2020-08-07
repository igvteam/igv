setlocal
::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

if exist %BatchPath%\jdk-11 (
  echo "Using bundled JDK."
  set JAVA_HOME=%BatchPath%\jdk-11
  set JAVA_CMD=%BatchPath%\jdk-11\bin\java
) else (
  echo "Using system JDK."
  set JAVA_CMD=java
)

if exist "%USERPROFILE%\.igv\java_arguments" (
  start %JAVA_CMD% -showversion --module-path=%BatchPath%\lib -Xmx1500m @%BatchPath%\igv.args @"%USERPROFILE%\.igv\java_arguments" --module=org.igv/org.broad.igv.tools.IgvTools %*
) else (
  start %JAVA_CMD% -showversion --module-path=%BatchPath%\lib -Xmx1500m @%BatchPath%\igv.args @"%USERPROFILE%\.igv\java_arguments" --module=org.igv/org.broad.igv.tools.IgvTools %*
)
