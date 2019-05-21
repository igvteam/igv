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

start %JAVA_CMD% -showversion --module-path=%BatchPath%\lib -Xmx1500m @igv.args --module=org.igv/org.broad.igv.tools.IgvTools gui
