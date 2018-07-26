=======================
IGV BINARY DISTRIBUTION
=======================

Prerequisites:

Java 9 or 10 (http://www.java.com).  Not compatible with Java 8 or Java 11 EA.


Instructions:

1. Download and unzip the distribution file to a directory of your choice.

2. To start IGV execute the following from the command line,

     java --module-path=lib -Xmx4g @igv.args --module org.broad.igv/org.broad.igv.ui.Main

Alternatively, you can start IGV with one of the following scripts.  Some of these may not
be present depending on the distribution you downloaded.  You might have to make the script 
executable (chmod a+x igv.sh).  


igv.bat       (for Windows)
igv.sh        (for LINUX and MAC OsX)
igv.command   (for MAC OsX, double-click to start)

The bat and shell scripts are configured to start IGV with 4GB of
memory.  This is a reasonable default for most machines.  If you are
working with very large datasets you can increase the amount of memory
available to IGV by editing the first line of the startup script.
Specifically change the value of the "-Xmx" parameter.  For example,
to start IGV with 8 gigabyte of memory  change the value

   -Xmx4g

to

   -Xmx8g

