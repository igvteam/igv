#!/bin/sh
cd `dirname $0`
java -Xmx1500m -Dapple.laf.useScreenMenuBar=true -jar `dirname $0`/lib/igvtools.jar gui
