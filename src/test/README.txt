
========================
TESTING IGV
========================

Prerequisities:
In addition to those used for building IGV,
additional data must be downloaded.


These are available for download at ftp://ftp.broadinstitute.org/pub/igv/largedata.zip.
This zip file is approximately 300 Mb.

Unzip these files to igv/test/largedata. The unzipped size is approximately 1 Gb.

If you push to put these files at a different path, this path must be passed to ant as the property LARGE_DATA_DIR.
For example, let's say you placed these data files in /user/remote/drive/largedata:

./gradlew -DLARGE_DATA_DIR="/user/remote/drive/largedata/" build

Also, some tests require a local instance of the Mongo executable.  You should pass this to the build as the
property MONGO_EXEC_PATH, for example:
./gradlew -DLARGE_DATA_DIR="/user/remote/drive/largedata/" -DMONGO_EXEC_PATH="/path/to/mongodb-2.4.6/bin/mongod" build