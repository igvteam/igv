setlocal

start jdk-11\bin\javaw -showversion --module-path=lib -Xmx8g -Dproduction=true @igv.args -Djava.net.preferIPv4Stack=true -Dsun.java2d.noddraw=true --module=org.igv/org.broad.igv.ui.Main 
