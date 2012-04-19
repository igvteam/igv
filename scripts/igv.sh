#apple.laf.useScreenMenuBar for Macs, to put menu bar at top of screen
#-Xmx750m indicates 750 mb of memory, adjust number up or down as desired
#Script must be in the same directory as igv.jar
java -Dapple.laf.useScreenMenuBar=true -Xmx750m -jar `dirname $0`/igv.jar $*