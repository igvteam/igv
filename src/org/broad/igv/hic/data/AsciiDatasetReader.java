package org.broad.igv.hic.data;

import org.broad.igv.hic.AlignmentsParser;

import java.io.*;

/**
 * @author jrobinso
 * @date Aug 18, 2010
 */
public class AsciiDatasetReader implements DatasetReader {

    String genomeId;
    File file;

    public AsciiDatasetReader(String genomeId, File f) {
        this.file = f;
        this.genomeId = genomeId;

    }

    public Dataset read() throws IOException {

        Dataset dataset = new Dataset(this);
        dataset.setCaching(false);
        return dataset;
    }

    public Matrix readMatrix(String key) throws IOException {

        String [] chroms = key.split("_");
        if(chroms.length != 2) {
            // todo error
            System.err.println("Uninterpretable matrix key: " + key);
            return null;
        }
        else {
            Integer chr1 = Integer.parseInt(chroms[0]);
            Integer chr2 = Integer.parseInt(chroms[1]);

            InputStream is = null;
            Matrix m;
            try {
                is = new FileInputStream(file);
                m = AlignmentsParser.readMatrix(is, genomeId, chr1, chr2);
            } finally {
                if(is != null) is.close();
            }

            return m;
        }

    }

}
