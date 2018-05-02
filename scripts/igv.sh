#!/bin/sh

#This script is intended for launch on *nix machines

#-Xmx4000m indicates 4000 mb of memory, adjust number up or down as needed
#Script must be in the same directory as igv.jar
#Add the flag -Ddevelopment = true to use features still in development
#Add the flag -Dsun.java2d.uiScale=2.0 for HiDPI displays
prefix=`dirname $(readlink $0 || echo $0)`
exec java --module-path="$prefix"/lib -Xmx4000m \
    @igv.args \
	-Dapple.laf.useScreenMenuBar=true \
	-Djava.net.preferIPv4Stack=true \
	--module org.broad.igv/org.broad.igv.ui.Main "$@"
