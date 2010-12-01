#!/bin/sh
cd `dirname $0`
java -Xmx1500m  -jar `dirname $0`/igvtools.jar --gui
