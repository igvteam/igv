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

The bat and shell scripts are configured to start IGV with 4GB of memory.  This is a 
reasonable default for most machines but if you are working with very large datasets
you can override this setting (and other Java-related defaults) by editing IGV's
java_arguments file, found here (create it if it doesn't exist):
   $HOME/.igv/java_arguments           (Mac and Linux)
   %USERPROFILE%/.igv/java_arguments   (Windows)

Specifically set the value of the "-Xmx" parameter, which will be commented-out (with
a '#' character) by default when this file is created.  For example, to start IGV with 
8 GB of memory, uncomment this line and set the value to 

   -Xmx8g

This will override the default 4GB memory specification.

Other Java-related command-line options can also be set in this file, though changing anything
beyond the memory specification is for advanced users only and is not recommended.  See
   https://docs.oracle.com/en/java/javase/11/tools/java.html
for more information on the Java 11 command line, and
   https://docs.oracle.com/en/java/javase/11/tools/java.html#GUID-4856361B-8BFD-4964-AE84-121F5F6CF111
in particular for specifics of the "java_arguments" file format.

The igv_hidpi.sh script is set up for 2x scaling.  To modify it to do 4x scaling, for 
example, change the value

   -Dsun.java2d.uiScale=2

to

   -Dsun.java2d.uiScale=4


Fractional values are *NOT* supported at this time.  Note that here again, you can add
this specification to the java_arguments file instead of editing the launcher scripts.
Doing so will allow the standard 'igv.sh' to work properly on HiDPI screens without
the need for the 'igv_hidpi.sh' script.
