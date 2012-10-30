#apple.laf.useScreenMenuBar for Macs, to put menu bar at top of screen
#-Xmx1200m indicates 1200 mb of memory, adjust number up or down as desired
#Script must be in the same directory as igv.jar
java -Dapple.laf.useScreenMenuBar=true -Xmx1200m -Djava.net.preferIPv4Stack=true -jar `dirname $0`/igv.jar $*