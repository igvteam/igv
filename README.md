# igv
[![Build Status](https://travis-ci.org/igvteam/igv.svg?branch=master)](https://travis-ci.org/igvteam/igv)

Integrative Genomics Viewer - desktop genome visualization tool for Mac, Windows, and Linux.

### Building

These instructions are meant for developers interested in working on the IGV code.  For normal use,
we recommend the pre-built releases available at [http://software.broadinstitute.org/software/igv/download](http://software.broadinstitute.org/software/igv/download).

Builds are executed from the IGV project directory.  Files will be created in the 'build' subdirectory.

IGV has been tested on **Java 11**. Previous (versions =< 2.6.3) running on Java8 have been deprecated.

NOTE: If on a Windows platform use ```./gradlew.bat'``` in the instructions below

#### Folder structure and build targets

Both [OpenJDK](https://openjdk.java.net/) and [Amazon's Correto Java 11](https://aws.amazon.com/corretto/) distributions have been tested with IGV.

* Install Gradle for your platform.  See https://gradle.org/ for details.

* Use ```./gradlew createDist``` to build a distribution directory (found in ```build/distributions```) containing 
  the igv.jar and its required runtime third-party dependencies (batik-codec, goby, and log4j-core) as
  well as helper scripts for launching.

    * These four JARs will be identical to those available in the download bundles from our website, 
    with the exception that they will not be signed with our certificate (required for JNLP) and
    will have slightly different build properties (timestamp, etc) in about.properties.

    * All four JARs must be in the same location in order to run IGV.  It can be run directly from
    'build/IGV-<YOUR_PLATFORM>'.

    * Launch IGV with `igv.sh` on Linux, `igv.command` on Mac, and `igv.bat`.
     These scripts can be edited to adjust JVM flags like maximum memory, etc.

    * To run igvtools from the command line use the script `igvtools` on Linux and Mac, or igvtools.bat
    on Windows. 

  
* Use ```./gradlew test``` to run the test suite.  See 'src/test/README.txt' for more information about running
  the tests.

Note that Gradle creates a number of other subdirectories in 'build'.  These can be safely ignored.

#### Amazon Web Services support

For more details on how to use IGV with AWS, please refer to the [UMCCR documentation](https://umccr.org/blog/igv-amazon-backend-setup/).
