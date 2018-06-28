::Get the current batch file's short path

for %%x in (%0) do set BatchPath=%%~dpsx

for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

java --module-path=%BatchPath%\lib -Xmx4000m -Dproduction=true @%BatchPath%\igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true --module org.broad.igv/org.broad.igv.ui.Main  %*
