The IGV JAR build now uses Gradle instead of Ant.  

Install Gradle for your platform.  See https://gradle.org/ for details.  At this time, *only* Java 8
is supported; we're working on support for Java 9+ in a future update.

Builds are executed from the IGV project directory.  Files will be created in the 'build' subdirectory.
You may need to execute 'gradle wrapper' to set up the gradle wrapper.  This should be necessary the 
first time only, or if you clean up the local .gradle directory.  Only do this if the './gradlew' calls
fail.

There are other options but these cover the most common uses:
- Use './gradlew createDist' to build a "distribution" directory (found in 'build/dist') containing 
  igv.jar, all of the third-party library JARs as well as helper scripts for launching.  
  Launch with 'igv.sh' on UNIX, 'igv.command' on Mac, and 'igv.bat' on Windows.  These scripts can be
  edited to adjust e.g. JVM flags like maximum memory, etc.
  Note that Gradle creates a number of other subdirectories in 'build'.  These can be safely ignored.
- Use './gradlew build' to build everything and run the test suite.  See 'src/test/README.txt' for more
  information about running the tests.
- OPTIONAL: assuming Ant is installed and configured with BCEL, a reduced-footprint version of the JAR
  can be produced by running './gradlew jar' or './gradlew build' and then 'ant -f build-shrink.xml'.  
  This JAR will be found in 'build/IGV-dist' along with two supporting JARs (batik-codec and goby) and
  helper scripts to run them (as above).
  All three JARs must be in the same location in order to run IGV. 
  These three JARs will be identical to those available in the download bundles from our website, with 
  the exception that they will not be signed with our certificate (required for JNLP).  


NOTE: In the above, use './gradlew.bat' on the Windows platform.