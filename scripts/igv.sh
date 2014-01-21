#!/bin/sh

#This script is intended for launch on *nix machines

#-Xmx2000m indicates 2000 mb of memory, adjust number up or down as needed
#-Dproduction=true disables non-released and development features
#Script must be in the same directory as igv.jar
#Add the flag -Ddevelopment = true to use features still in development
prefix=`dirname $(readlink $0 || echo $0)`
exec java -Xmx2000m \
	-Dapple.laf.useScreenMenuBar=true \
	-Djava.net.preferIPv4Stack=true \
	-jar "$prefix"/igv.jar "$@"
