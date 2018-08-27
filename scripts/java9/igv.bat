::Get the current batch file's short path

for %%x in (%0) do set BatchPath=%%~dpsx

for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

java --module-path=%BatchPath%\lib -Xmx4g -Dproduction=true @%BatchPath%\igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true --module org.igv/org.broad.igv.ui.Main  %*
