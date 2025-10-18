#!/bin/sh

#This script is intended for launch on *nix machines

#-Xmx8g indicates 8 gb of memory.
#To adjust this (or other Java options), edit the "$HOME/.igv/java_arguments" 
#file.  For more info, see the README at 
#https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt 
#Add the flag -Ddevelopment = true to use features still in development
#Add the flag -Dsun.java2d.uiScale=2 for HiDPI displays
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
        -Dsamjdk.snappy.disable=true \
        -Dapple.laf.useScreenMenuBar=true \
        -Djava.net.preferIPv4Stack=true \
        -Djava.net.useSystemProxies=true \
        @"$HOME/.igv/java_arguments" \
        -cp "$CP" \
        org.broad.igv.ui.Main "$@"
else
    java -Xmx8g \
        @"${prefix}/igv.args" \
        -Dsamjdk.snappy.disable=true \
        -Dapple.laf.useScreenMenuBar=true \
        -Djava.net.preferIPv4Stack=true \
        -Djava.net.useSystemProxies=true \
        -cp "$CP" \
        org.broad.igv.ui.Main "$@"
fi