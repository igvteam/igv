#!/bin/sh
#apple.laf.useScreenMenuBar for Macs, to put menu bar at top of screen
#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed
#-Dproduction=true disables non-released and development features
#Script must be in the same directory as igv.jar
prefix=`dirname $(readlink $0 || echo $0)`
exec java -Xmx2000m \
	-Dproduction=true \
	-Dapple.laf.useScreenMenuBar=true \
	-Djava.net.preferIPv4Stack=true \
	-jar "$prefix"/igv.jar "$@" &
