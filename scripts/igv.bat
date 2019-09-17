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

start %JAVA_CMD% -showversion --module-path=%BatchPath%\lib -Xmx4g -Dproduction=true @%BatchPath%\igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true --module=org.igv/org.broad.igv.ui.Main  %*
