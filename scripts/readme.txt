=======================
IGV BINARY DISTRIBUTION
=======================

Prerequisites:

Java 8.0 (http://www.java.com).  Not compatible with Java 9+


Instructions:

1. Download and unzip the distribution file to a directory of your choice.

2. To start IGV execute the following from the command line,

     java -Xmx750m -jar igv.jar

Alternatively, you can start IGV with one of the following scripts.  You might have to make the script executable
(chmod a+x igv.sh).


igv.bat       (for Windows)
igv.sh        (for Linux and macOS)
igv.command   (for macOS, double-click to start)

The shell scripts are configured to start IGV with 1500MB of memory (1 GB
for the bat script).  This is a reasonable default for most machines.  If 
you are working with very large datasets you can increase the amount of 
memory available to IGV by editing the first line of the startup script.
Specifically change the value of the "-Xmx" parameter.  For example,
to start IGV with 1 GB of memory  change the value

   -Xmx1500m

to

   -Xmx1000m

