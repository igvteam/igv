::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx
java -Xmx750m -Dsun.java2d.noddraw=true -jar %BatchPath%\igv.jar  %*
