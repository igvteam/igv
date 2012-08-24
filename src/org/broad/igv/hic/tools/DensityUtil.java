package org.broad.igv.hic.tools;

import org.apache.log4j.Logger;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.HiCGlobals;
import org.broad.igv.hic.data.DensityFunction;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 1/1/12
 */
public class DensityUtil {
    private static Logger log = Logger.getLogger(DensityCalculation.class);

    public static void main(String[] args) throws IOException {

        String genomeID = "b37";
        List<Chromosome> chromosomes = HiCTools.loadChromosomes(genomeID);
        calculate(chromosomes);
        //read("/xchip/igv/dev/hic/testFiles/test.hic.densities");
        // dumpDensities("/xchip/igv/dev/hic/testFiles/test.hic.densities", 1, 14); //Hi-C_HindIII_Human_August.hic.densities", 1, 14);
    }

    /*private static void dumpDensities(String path, int zoomNumber, int chr) throws IOException {
        InputStream is = ParsingUtils.openInputStream(path);
        Map<Integer, DensityFunction> d = readDensities(is);
        DensityFunction df = d.get(zoomNumber);
        if (df != null) {
            for (int b = 0; b < 250; b++) {
                System.out.println(df.getDensity(chr, b));
            }
        }
        is.close();

    }

    private static Map<Integer, DensityFunction> read(String ifile) throws IOException {
        InputStream is = ParsingUtils.openInputStream(ifile);
        return readDensities(is);
    }
*/
    private static void calculate(List<Chromosome> chromosomes) throws IOException {

        String[] paths = {"/xchip/igv/dev/hic/testFiles/GSM455139_428EGAAXX.7.maq.hic.summary.binned.txt", "/xchip/igv/dev/hic/testFiles/GSM455140_428EGAAXX.8.maq.hic.summary.binned.txt"};


        Map<String, Integer> chrIndexMap = new HashMap<String, Integer>();
        for (Chromosome chr : chromosomes) {
            if (chr != null && chr.getIndex() > 0)
                chrIndexMap.put(chr.getName(), chr.getIndex());
        }


        // Limit calcs to 10KB
        int[] gridSizeArray = new int[8];
        for (int i = 0; i < 8; i++) {
            gridSizeArray[i] = HiCGlobals.zoomBinSizes[i];
        }

        DensityCalculation[] calcs = new DensityCalculation[gridSizeArray.length];
        for (int z = 0; z < gridSizeArray.length; z++) {
            calcs[z] = new DensityCalculation(chromosomes, gridSizeArray[z]);
        }

        for (String path : paths) {
            PairIterator iter = (path.endsWith(".bin")) ?
                    new BinPairIterator(path) :
                    new AsciiPairIterator(path, chrIndexMap);
            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                if (pair.getChr1() == (pair.getChr2())) {

                    int index = pair.getChr1();
                    for (int z = 0; z < gridSizeArray.length; z++) {
                        calcs[z].addDistance(index, pair.getPos1(), pair.getPos2());
                    }
                }

            }
        }
        for (int z = 0; z < gridSizeArray.length; z++) {
            calcs[z].computeDensity();
        }

        //outputDensities(calcs, new File("/xchip/igv/dev/hic/testFiles/HindIII_Human_August.densities"));
    }


/*    private static void outputDensities(DensityCalculation[] calcs, File outputFile) throws IOException {

        LittleEndianOutputStream os = null;
        try {
            os = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            os.writeInt(calcs.length);
            for (int i = 0; i < calcs.length; i++) {
                calcs[i].outputBinary(os);
            }
        } finally {
            if (os != null) os.close();
        }
    }*/

    /**
     * Return a map of zoom level -> DensityFunction
     *
     * @param les
     * @return
     * @throws IOException
     */
    public static Map<Integer, DensityFunction> readDensities(LittleEndianInputStream les, boolean isNewVersion) throws IOException {

        int nZooms = les.readInt();
        Map<Integer, DensityFunction> densityMap = new HashMap<Integer, DensityFunction>();
        // TODO -- Its assumed densities are in number order and indeces match resolutions.  This is fragile,
        // encode resolutions in the next round
        for (int i = 0; i < nZooms; i++) {
            DensityCalculation calc = new DensityCalculation(les, isNewVersion);
            densityMap.put(i, new DensityFunction(calc));
        }

        return densityMap;

    }

    public static Map<Integer, DensityFunction> readDensities(LittleEndianInputStream les) throws IOException {
        return readDensities(les, true);
    }
}