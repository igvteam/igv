package org.broad.igv.hic.data;

import org.broad.igv.hic.tools.Preprocessor;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * @author jrobinso
 *         Date: 10/17/12
 *         Time: 8:38 AM
 */
public interface DatasetReader {
    String getPath();

    Dataset read() throws FileNotFoundException;

    int getVersion();

    Matrix readMatrix(String key) throws IOException;

    Block readBlock(int blockNumber, Preprocessor.IndexEntry idx) throws IOException;
}
