package org.broad.igv.hic.data;

import org.broad.igv.feature.Chromosome;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 12, 2010
 */
public class Dataset {

    // this needs to be read in from the file instead of left as a global
    // will need to be part of the header.
    private static final int[] zoomBinSizes = {2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1}; //,  2500, 1000};
    private static final String[] zoomLabels = {"2.5 MB", "1 MB", "500 KB", "250 KB", "100 KB", "50 KB", "25 KB", "10 KB", "5 KB", "frag"};//, "2.5 KB", "1 KB"};
    private static final int nZooms = 10;

    private boolean caching = true;

    //Chromosome lookup table
    public  Chromosome [] chromosomes;

    Map<String, Matrix> matrices = new HashMap<String, Matrix>(25 * 25);

    private DatasetReader reader;
    private Map<Integer, DensityFunction> df;

    public Dataset(DatasetReader reader) {
        this.reader = reader;
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

    public static int getNumberZooms() {
        return nZooms;
    }

//    public static int getNumberZoomBinsToWrite(){
//        return zoomLabels.length;
//    }


    public static String getZoomLabel(int index) {
        return zoomLabels[index];
    }

    public static int getZoom(int index) {
        return zoomBinSizes[index];
    }

    public Map<Integer, DensityFunction> getZoomToDensity() {
        return df;
    }

    public void setZoomToDensity(Map<Integer, DensityFunction> df) {
        this.df = df;
        for (Map.Entry<Integer, DensityFunction> entry : df.entrySet())
        {
            entry.getValue().setChromosomes(chromosomes);
        }
    }

    public Chromosome[] getChromosomes() {
        return chromosomes;
    }

    public void setChromosomes(Chromosome[] chromosomes) {
        this.chromosomes = chromosomes;
    }

    public int getVersion() {
        return reader.getVersion();
    }
}
