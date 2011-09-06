package org.broad.igv.hic.data;

import java.io.File;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Aug 18, 2010
 */
public interface DatasetReader {

    Dataset read() throws IOException;
    
    Matrix readMatrix(String key) throws IOException;

}
