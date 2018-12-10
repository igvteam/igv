# igv
[![Build Status](https://travis-ci.org/igvteam/igv.svg?branch=master)](https://travis-ci.org/igvteam/igv)

Integrative Genomics Viewer - dekstop genome visualization tool for Mac, Windows, and Linux.

### Building

These instructions are meant for developers interested in working on the IGV code.  For normal use,
we recommend the pre-built releases available at [http://software.broadinstitute.org/software/igv/download](http://software.broadinstitute.org/software/igv/download).


Builds are executed from the IGV project directory.  Files will be created in the 'build' subdirectory.
You may need to execute 'gradle wrapper' to set up the gradle wrapper.  This should be necessary the 
first time only, or if you clean up the local .gradle directory.  Only do this if the './gradlew' calls
fail.

There are two different build options, one for Java 8 and another for Java 9 and up.  The default is 
to build for Java 8.  Java 8 builds are *NOT* compatible with Java 9 and vice versa.  

There are other options but these cover the most common uses:

NOTE: If on a Windows platform use ```./gradlew.bat'``` in the instructions below

#### Java 8

* Install Gradle for your platform.  See https://gradle.org/ for details.

* Use ```./gradlew createDist``` to build a distribution directory (found in ```build/IGV-dist```) containing 
  the igv.jar and its required runtime third-party dependencies (batik-codec, goby, and log4j-core) as
  well as helper scripts for launching.
  
    * These four JARs will be identical to those available in the download bundles from our website, 
    with the exception that they will not be signed with our certificate (required for JNLP) and
    will have slightly different build properties (timestamp, etc) in about.properties.
    * All four JARs must be in the same location in order to run IGV.  It can be run directly from
    'build/IGV-dist' 
    
  *  Launch with 'igv.sh' on UNIX, 'igv.command' on Mac, and 'igv.bat' on Windows.  These scripts can
    be edited to adjust JVM flags like maximum memory, etc.
    
  *  All other runtime dependencies are bundled into igv.jar.  There is also an igv-minimal.jar in
    'build/libs' containing just the IGV classes and resources for those who prefer to manage 
    dependencies as separate files.

    
* Use ```./gradlew createToolsDist``` to build an igvtools distribution directory (found in 
  'build/IGVTools-dist') containing the igvtools.jar and dependencies (same as for IGV, above)  
  JAR dependencies plus helper scripts for running and launching.
  As above, these JARs will be identical aside from signing, timestamps, etc. and all must be
  present together to run.  See igvtools_readme.txt for more info.
  
* Use ```./gradlew test``` to run the test suite.  See 'src/test/README.txt' for more information about running
  the tests.

Note that Gradle creates a number of other subdirectories in 'build'.  These can be safely ignored.

#### Java 9

The instructions for Java 9 are nearly identical other than the need to specify the Java 9 build file
and that the results will be found in 'build_java9' rather than 'build'.  More specifically:

* Use './gradlew -b build_java9.gradle createDist' to build a distribution directory with helper scripts
  for launching.  The structure is slightly different but the concept is the same.
  
* Use './gradlew -b build_java9.gradle createToolsDist' for the igvtools distribution.

* Use './gradlew -b build_java9.gradle test' to run the test suite.

The full JAR build option is *NOT* available for Java 9+ because of modularity requirements.

