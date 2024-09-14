#!/bin/sh
cd `dirname $0`
prefix=`dirname $(readlink -f $0 || echo $0)`

# Check whether or not to use the bundled JDK
if [ -d "${prefix}/jdk-21" ]; then
    echo echo "Using bundled JDK."
    JAVA_HOME="${prefix}/jdk-21"
    PATH=$JAVA_HOME/bin:$PATH
else
    echo "Using system JDK."
fi

# Check if there is a user-specified Java arguments file
if [ -e "$HOME/.igv/java_arguments" ]; then
    java -showversion --module-path="${prefix}/lib" -Xmx1500m \
        @"${prefix}/igv.args" \
        -Dapple.laf.useScreenMenuBar=true \
        @"$HOME/.igv/java_arguments" \
        --module=org.igv/org.broad.igv.tools.IgvTools gui
else
    java -showversion --module-path="${prefix}/lib" -Xmx1500m \
        @"${prefix}/igv.args" \
        -Dapple.laf.useScreenMenuBar=true \
        --module=org.igv/org.broad.igv.tools.IgvTools gui
fi
