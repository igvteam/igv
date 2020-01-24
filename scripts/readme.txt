=======================
IGV BINARY DISTRIBUTION
=======================

Prerequisites:

Java 11 (http://openjdk.java.net).  This is bundled with our distributions.
Not compatible with Java 8, 9, 10.


Instructions:

1. Download and unzip the distribution file to a directory of your choice.

2. To start IGV execute the following from the command line,

     java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main

Note that the command line has become more complex with Java 11 compared to Java 8.  
Alternatively, you can start IGV with one of the following scripts; this is 
recommended.  Some of these may not be present depending on the distribution you 
downloaded.  You might have to make the script executable (chmod a+x igv.sh).  


igv-launcher.bat  (for Windows)
igv.bat           (for Windows batch jobs)
igv.sh            (for Linux and macOS)
igv_hidpi.sh      (for Linux with HiDPI displays)
igv.command       (for macOS, double-click to start)

The bat and shell scripts are configured to start IGV with 4GB of
memory.  This is a reasonable default for most machines.  If you are
working with very large datasets you can increase the amount of memory
available to IGV by editing the first line of the startup script.
Specifically change the value of the "-Xmx" parameter.  For example,
to start IGV with 8 GB of memory change the value

   -Xmx4g

to

   -Xmx8g

The igv_hidpi.sh script is set up for 2x scaling.  To modify it to do 4x scaling, for 
example, change the value

   -Dsun.java2d.uiScale=2

to

   -Dsun.java2d.uiScale=4

Fractional values are *NOT* supported at this time. 