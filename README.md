# igv
[![Build Status](https://travis-ci.com/igvteam/igv.svg?branch=master)](https://travis-ci.com/igvteam/igv)
![GitHub issues](https://img.shields.io/github/issues/igvteam/igv)
![GitHub closed issues](https://img.shields.io/github/issues-closed/igvteam/igv)
![](https://img.shields.io/npm/l/igv.svg)

Integrative Genomics Viewer - desktop genome visualization tool for Mac, Windows, and Linux.

### Building

These instructions are meant for developers interested in working on the IGV code.  For normal use,
we recommend the pre-built releases available at [http://software.broadinstitute.org/software/igv/download](http://software.broadinstitute.org/software/igv/download).

Builds are executed from the IGV project directory.  Files will be created in the 'build' subdirectory.

IGV is tested on **Java 11**. Previous (versions =< 2.6.3) running on Java8 have been deprecated.

NOTE: If on a Windows platform use ```./gradlew.bat'``` in the instructions below

#### Folder structure and build targets

The IGV bundles ship with embedded JREs from AdoptOpenJDK.

* Install Gradle for your platform.  See https://gradle.org/ for details.

* Use ```./gradlew createDist``` to build a distribution directory (found in ```build/IGV-dist```) containing 
  the igv.jar and its required runtime third-party dependencies as well as helper scripts for launching.

    * Launch IGV with `igv.sh` or `igv_hidpi.sh` on Linux, `igv.command` on Mac, and `igv.bat` on Windows.

    * To run igvtools from the command line use the script `igvtools` on Linux and Mac, or igvtools.bat
      on Windows.  See the instructions in igvtools_readme.txt in that directory.

    * The launcher scripts expect this folder structure in order to run IGV.
  
* Use ```./gradlew test``` to run the test suite.  See 'src/test/README.txt' for more information about running
  the tests.
  
* This dashboard describes [project structure and dependencies](https://sourcespy.com/github/igvteamigv/). New contributors can quickly grasp overall structure of the code and technologies involved.  

* See this [README](https://raw.githubusercontent.com/igvteam/igv/master/scripts/readme.txt) for tips about using the IGV launcher scripts.

Note that Gradle creates a number of other subdirectories in 'build'.  These can be safely ignored.

#### Amazon Web Services support

Public data files hosted in Amazon S3 buckets can be loaded into IGV using [https endpoints](https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingBucket.html).  Authenticated access can be configured using the UMCCR contributed AWS configuration option.  For more details, refer to the [UMCCR documentation](https://umccr.org/blog/igv-amazon-backend-setup/).   Note authenticated access is currently restricted to the following file formats:  BAM, CRAM, tabix indexed VCF (vcf.gz / vcf.gz.tbi), and non-indexed feature file formats.
