
This archive contains the source distribution for the Integrative Genomics
Viewer (IGV).   


========================
BUILDING IGV
========================


Prerequisites: 

Java J2SE 6.0 or greater (http://java.sun.com/javase/download)
Ant 1.7.0 or greater (http://ant.apache.org/)
Optional: BCEL (http://commons.apache.org/bcel/)


1.  Download and unzip the source distribution file.  

2.  Run the provided ant script by running "ant" from the root directory
    of the distribution.


The above script will build "igv.jar" in the root directory of the distribution.
This file does not include the external libraries on which IGV depends, and so will only when located
in the same directory as a folder called "lib" containing all necessary jars.

Optional:

To package IGV for distribution, run the following command:

ant -Dinclude.libs=true

The command

ant -Dinclude.libs=false

will have the same effect as running "ant".


This requires the BCEL library.
Some ant distributions do not come with the optional library needed
needed to use BCEL. BCEL is used only in certain build tasks, and is
not needed for the default build. It is used when publishing the jar,
to remove unnecessary files for minimal file size.

When this library is absent, you may get an error
message like this:
Unable to load dependency analyzer:
org.apache.tools.ant.util.depend.bcel.FullAnalyzer

To solve this issue, download the latest version of ant from Apache, and make sure
it includes ant-apache-bcel.jar and bcel-5.2.jar (or latest version number).

========================
TESTING IGV
========================

Prerequisites:
In addition to those needed for building IGV, additional test
data must be downloaded. See test/README.txt for details and instructions.

===========================
RUNNING
===========================

After building igv.jar IGV can be launched by executing one of the following
command line scripts:

igv.bat       (for Windows systems)
igv.sh        (for LINUX and MAC OsX)
igv.command    (for MAC OsX, double-click to start)

The bat and shell scripts are configured to start IGV with 750MB of
memory.  This is a reasonable default for most machines.  If you are
working with very large datasets you can increase the amount of memory
available to IGV by editing the first line of the startup script.
Specifically change the value of the "-Xmx" parameter.  For example,
to start IGV with 1 gigabyte of memory  change the value

   -Xmx750m

to

   -Xmx1000m


