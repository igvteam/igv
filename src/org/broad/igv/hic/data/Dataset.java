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

    //Chromosome lookup table
    public static Chromosome [] chromosomes = new Chromosome [] {
        new Chromosome(0, "2L", 23011544),
        new Chromosome(1, "2R", 21146708),
        new Chromosome(2, "3L", 24543557),
        new Chromosome(3, "3R", 27905053),
        new Chromosome(4, "4", 1351857),
        new Chromosome(5, "X", 22422827),
        new Chromosome(6, "U", 10049037)
    };

    Map<String, Matrix> matrices = new HashMap(25 * 25);

    private DatasetReader reader;

    public Dataset(DatasetReader reader) {
        this.reader = reader;
        readChromosomes();
    }

    public void readChromosomes() {
        // todo -- read the chromosome dictionary for the fil
    }

    public Matrix getMatrix(Chromosome chr1, Chromosome chr2) {

        // order is arbitrary, convention is lower # chr first
        int t1 = Math.min(chr1.getIndex(), chr2.getIndex());
        int t2 = Math.max(chr1.getIndex(), chr2.getIndex());

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

    public Chromosome[] getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(Chromosome[] chromosomes) {
        this.chromosomes = chromosomes;
    }
}
