::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx
java --module-path=%BatchPath%\lib -Xmx1500m @igv.args --module org.igv/org.broad.igv.tools.IgvTools gui
