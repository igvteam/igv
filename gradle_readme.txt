The IGV JAR build now uses Gradle instead of Ant.  

Install Gradle for your platform.  See https://gradle.org/ for details.

Builds are executed from the IGV project directory.  Files will be created in the 'build' subdirectory.
You may need to execute 'gradle wrapper' to set up the gradle wrapper.  This should be necessary the 
first time only, or if you clean up the local .gradle directory.  Only do this if the './gradlew' calls
fail.

There are other options but these cover the most common uses:
- Use './gradlew jar' to build a runnable JAR, which will be created in 'build/lib'.  To run, this JAR 
  expects the third-party library JARs to be present in the same directory.
- Use './gradlew distZip' to build a distribution ZIP containing the runnable JAR plus all of the
  third-party library JARs and helper scripts for launching.
- Use './gradlew build' to build everything and run the test suite.  See 'src/test/README.txt' for more
  information about running the tests.

NOTE: In the above, use './gradlew.bat' if you're on Windows.