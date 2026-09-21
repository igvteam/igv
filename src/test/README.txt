
========================
TESTING IGV
========================

Prerequisities:
Those used for building IGV.  Test data is in igv/test/data.

Some tests require a local instance of the Mongo executable.  You should pass this to the build as the
property MONGO_EXEC_PATH, for example:
./gradlew -DMONGO_EXEC_PATH="/path/to/mongodb-2.4.6/bin/mongod" build
