package org.broad.igv.hic.data;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.HiC;
import org.broad.igv.hic.tools.Preprocessor;

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
    public Chromosome[] chromosomes;

    Map<String, Matrix> matrices = new HashMap<String, Matrix>(25 * 25);

    private DatasetReader reader;
    private Map<String, DensityFunction> df;
    private String genomeId;

    private int[] bpBinSizes;
    private int[] fragBinSizes;


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

                if (caching) matrices.put(key, m);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        return m;

    }

    public int getNumberZooms(HiC.Unit unit) {
        return  unit == HiC.Unit.BP ?  bpBinSizes.length : fragBinSizes.length;
    }

    public int getZoom(HiC.Unit unit, int index) {
        return unit == HiC.Unit.BP ?  bpBinSizes[index] : fragBinSizes[index];
    }

    /**
     * Function needed for legacy datasets.
     *
     * @param zoom
     * @return
     */
    public DensityFunction getDensityFunction(int zoom) {

        if (df == null) return null;

        int binSize = Preprocessor.bpBinSizes[zoom];
        String unit = binSize == 1 ? "FRAG" : "BP";

        return getExpectedValues(unit, binSize);
    }

    public DensityFunction getExpectedValues(String unit, int binSize) {

        String key = unit + "_" + binSize;

        return df.get(key);
    }


    public void setDensityFunctionMap(Map<String, DensityFunction> df) {
        this.df = df;
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

    public void setGenomeId(String genomeId) {
        this.genomeId = genomeId;
    }

    public void setBpBinSizes(int[] bpBinSizes) {
        this.bpBinSizes = bpBinSizes;
    }

    public void setFragBinSizes(int[] fragBinSizes) {
        this.fragBinSizes = fragBinSizes;
    }

}
