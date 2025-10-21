#!/bin/sh

#This script is intended for launching on Macs
#It may or may not work on *nix, definitely not on windows

#-Xmx4g indicates 4 gb of memory.
#To adjust this (or other Java options), edit the "$HOME/.igv/java_arguments" 
#file.  For more info, see the README at 
#https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt 
#apple.laf.useScreenMenuBar for Macs, to put menu bar at top of screen
#-Xdock:name again for Macs, sets the name in menu bar
#Add the flag -Ddevelopment = true to use features still in development
prefix="$(cd "$(dirname "$0")" && pwd)"

# Check whether or not to use the bundled JDK
if [ -d "${prefix}/jdk-21" ]; then
    JAVA_HOME="${prefix}/jdk-21"
    PATH=$JAVA_HOME/bin:$PATH
else
    echo "Using system JDK. IGV requires Java 21." >&2
fi

# Report on Java version
java -version

# Build classpath including all jars in lib (non-recursive)
CP="${prefix}/lib/*"

# Check if there is a user-specified Java arguments file
if [ -e "$HOME/.igv/java_arguments" ]; then
  java -Xmx8g \
         @"${prefix}/igv.args" \
        -Xdock:name="IGV" \
        -Xdock:icon="${prefix}/IGV_64.png" \
        -Dapple.laf.useScreenMenuBar=true \
        -Djava.net.preferIPv4Stack=true \
        -Djava.net.useSystemProxies=true \
        -cp "$CP" \
        org.broad.igv.ui.Main "$@"
else
    java -showversion  -Xmx8g \
        @"${prefix}/igv.args" \
        -Xdock:name="IGV" \
        -Xdock:icon="${prefix}/IGV_64.png" \
        -Dapple.laf.useScreenMenuBar=true \
        -Djava.net.preferIPv4Stack=true \
        -Djava.net.useSystemProxies=true \
        -cp "$CP" \
        org.broad.igv.ui.Main "$@"
fi
