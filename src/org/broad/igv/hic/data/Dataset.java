package org.broad.igv.hic.data;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 12, 2010
 */
public class Dataset {

    private boolean caching = true;

    // TODO -- use soft reference cache
    Map<String, Matrix> matrices = new HashMap(25 * 25);

    private DatasetReader reader;

    public Dataset(DatasetReader reader) {
        this.reader = reader;
    }

    public Matrix getMatrix(int chr1, int chr2) {

        // order is arbitrary, convention is lower # chr first
        int t1 = Math.min(chr1, chr2);
        int t2 = Math.max(chr1, chr2);

        String key = Matrix.generateKey(t1, t2);
        Matrix m = matrices.get(key);

        if (m == null && reader != null) {
            try {
                m = reader.readMatrix(key);
                if(caching) matrices.put(key, m);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        return m;

    }


    public boolean isCaching() {
        return caching;
    }

    public void setCaching(boolean caching) {
        this.caching = caching;
    }
}
