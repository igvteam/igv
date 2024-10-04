#!/bin/bash

#This script is intended for launching on Macs
#It is intended to be compiled be 'shc' on Mac in order to have the
#launcher be a binary executable rather than a script.  You can obtain
#'sch' through MacPorts or Homebrew.  Then, the steps to compile look
#like this:
#  $ export CFLAGS=-mmacosx-version-min=10.10   # Pass-thru to 'cc' for backward compatibility.
#  $ shc -r -f IGV.sh  # The '-r' means "relaxed security" and is more portable.
#  $ mv IGV.sh.x IGV   # Rename the executable to the expected name.

#-Xmx8g indicates 8 gb of memory.
#To adjust this (or other Java options), edit the "$HOME/.igv/java_arguments" 
#file.  For more info, see the README at 
#https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt 
#Add the flag -Ddevelopment = true to use features still in development
prefix=`dirname $(readlink $0 || echo $0)`

# Check whether or not to use the bundled JDK
echo ${prefix}
if [ -d "${prefix}/../jdk-17" ]; then
    echo echo "Using bundled JDK."
    JAVA_HOME="${prefix}/../jdk-17"
    PATH=$JAVA_HOME/bin:$PATH
else
    echo "Using system JDK."
    java -version
fi

# Check if there is a user-specified Java arguments file
if [ -e "$HOME/.igv/java_arguments" ]; then
    EXTRA_ARGS_FILE=@"$HOME/.igv/java_arguments"
fi

exec java -showversion --module-path="${prefix}/../Java/lib" -Xmx8g \
    @"${prefix}/../Java/igv.args" \
    -Xdock:name="IGV" \
    -Xdock:icon="${prefix}/../Resources/IGV_64.png" \
    -Dapple.laf.useScreenMenuBar=true \
    -Djava.net.preferIPv4Stack=true ${EXTRA_ARGS_FILE} \
    --module=org.igv/org.broad.igv.ui.Main
