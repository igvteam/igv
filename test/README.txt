To run the tests, simply execute the command:
ant alltests

Some test data files are too large to be stored in the SVN repository.
These are available for download at ftp://ftp.broadinstitute.org/pub/igv/largedata.zip.
This zip file is approximately 300 Mb.

Unzip these files to igv/test/largedata. The unzipped size is approximately 1 Gb.

If you push to put these files at a different location, this path must be passed to ant as the property LARGE_DATA_DIR. For example, let's say you placed these data files in /user/remote/drive/largedata:

ant -DLARGE_DATA_DIR="/user/remote/drive/largedata" alltests