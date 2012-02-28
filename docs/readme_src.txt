
This archive contains the source distribution for the Integrative Genomics
Viewer (IGV).   


========================
BUILDING IGV
========================


Prerequisites: 

Java J2SE 6.0 or greater (http://java.sun.com/javase/download)
Ant 1.7.0 or greater (http://ant.apache.org/)


1.  Download and unzip the source distribution file.  

2.  Run the provided ant script by running "ant" from the root directory
    of the distribution.


The above script will build "igv.jar" in the root directory of the distribution.


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

