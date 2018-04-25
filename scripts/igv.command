#!/bin/sh

#This script is intended for launching on Macs
#It may or may not work on *nix, definitely not on windows

#apple.laf.useScreenMenuBar for Macs, to put menu bar at top of screen
#-Xdock:name again for Macs, sets the name in menu bar
#-Xmx4000m indicates 4000 mb of memory, adjust number up or down as needed
#Script must be in the same directory as igv.jar
#Add the flag -Ddevelopment = true to use features still in development
prefix=`dirname $(readlink $0 || echo $0)`
exec java  -classpath "$prefix"/lib --module-path="$prefix"/modules -Xmx4000m \
    -Xdock:name="IGV" \
	-Dapple.laf.useScreenMenuBar=true \
	-Djava.net.preferIPv4Stack=true \
	--module org.broad.igv/org.broad.igv.ui.Main "$@"
