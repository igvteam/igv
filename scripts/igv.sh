#!/bin/sh

#This script is intended for launch on *nix machines

#-Xmx4g indicates 4 gb of memory, adjust number up or down as needed
#Script must be in the same directory as igv.jar
#Add the flag -Ddevelopment = true to use features still in development
prefix=`dirname $(readlink $0 || echo $0)`
exec java -Xmx4g \
    -Dapple.laf.useScreenMenuBar=true \
    -Dawt.useSystemAAFontSettings=on -Dswing.aatext=true \
    -Djava.net.preferIPv4Stack=true \
    -jar "$prefix"/lib/igv.jar "$@"
