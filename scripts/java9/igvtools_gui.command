#!/bin/sh
cd `dirname $0`
java --module-path=`dirname $0`/lib -Xmx1500m \
    @igv.args \
    -Dapple.laf.useScreenMenuBar=true \
    --module org.broad.igv/org.broad.igv.tools.IgvTools gui
