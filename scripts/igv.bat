setlocal
::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

if exist %BatchPath%\jdk-21 (
  echo "Using bundled JDK."
  set JAVA_HOME=%BatchPath%\jdk-21
  set JAVA_CMD=%BatchPath%\jdk-21\bin\javaw
) else (
  echo "Using system JDK."
  set JAVA_CMD=java
)

::-Xmx8g indicates 8 gb of memory.
::To adjust this (or other Java options), edit the "%USERPROFILE%\.igv\java_arguments" 
::file.  For more info, see the README at 
::https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt 

set CP=%BatchPath%lib\*

if exist "%USERPROFILE%\.igv\java_arguments" (
    %JAVA_CMD% -Xmx8g ^
        @%BatchPath%igv.args ^
        -Dsamjdk.snappy.disable=true ^
        -Djava.net.preferIPv4Stack=true ^
        -Djava.net.useSystemProxies=true ^
        @"%USERPROFILE%\.igv\java_arguments" ^
        -cp "%CP%" org.broad.igv.ui.Main %*
) else (
    %JAVA_CMD% -Xmx8g ^
        @%BatchPath%igv.args ^
        -Dsamjdk.snappy.disable=true ^
        -Djava.net.preferIPv4Stack=true ^
        -Djava.net.useSystemProxies=true ^
        -cp "%CP%" org.broad.igv.ui.Main %*
)
