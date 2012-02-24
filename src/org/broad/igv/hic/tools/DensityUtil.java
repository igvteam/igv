package org.broad.igv.hic.tools;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broad.igv.hic.MainWindow;
import org.broad.igv.hic.data.Chromosome;
import org.broad.igv.hic.data.DensityFunction;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 1/1/12
 */
public class DensityUtil {
    private static Logger log = Logger.getLogger(DensityCalculation.class);
    
    public static void main(String[] args) throws IOException {
         calculate();
        //read();
        dumpDensities("/xchip/igv/dev/hic/testFiles/test.hic.densities", 1, 14); //Hi-C_HindIII_Human_August.hic.densities", 1, 14);
    }

    private static void dumpDensities(String path, int zoomNumber, int chr) throws IOException {
        InputStream is = ParsingUtils.openInputStream(path);
        Map<Integer, DensityFunction> d = readDensities(is);
        DensityFunction df = d.get(zoomNumber);
        for(int b=0; b<250; b++) {
            System.out.println(df.getDensity(chr, b));
        }


        is.close();

    }

    private static Map<Integer, DensityFunction> read(String ifile) throws IOException {
        //String ifile = "/Users/jrobinso/IGV/hic/test.hic.densitites";
        InputStream is = ParsingUtils.openInputStream(ifile);
        Map<Integer, DensityFunction> d = readDensities(is);
        return d;
    }

    private static void calculate() throws IOException {
//        String[] paths = {"/Users/jrobinso/IGV/hic/formattedalignment.txt.gz"};
        String[] paths = {"/xchip/igv/dev/hic/testFiles/GSM455139_428EGAAXX.7.maq.hic.summary.binned.txt"};
//                "/xchip/igv/dev/hic/testFiles/GSM455140_428EGAAXX.8.maq.hic.summary.binned.txt"};
//        //chromosomes = HiCTools.b37Chromosomes;
        Chromosome[] chromosomes = HiCTools.hg18Chromosomes;
//        Chromosome[] chromosomes = HiCTools.b37Chromosomes;
//        String[] paths = {
//                "/Volumes/igv/data/broad/hic/human/GSM455133_30E0LAAXX.1.maq.hic.summary.binned.txt.gz",
//                "/Volumes/igv/data/broad/hic/human/GSM455134_30E0LAAXX.2.maq.hic.summary.binned.txt.gz",
//                "/Volumes/igv/data/broad/hic/human/GSM455135_30U85AAXX.2.maq.hic.summary.binned.txt.gz",
//                "/Volumes/igv/data/broad/hic/human/GSM455136_30U85AAXX.3.maq.hic.summary.binned.txt.gz",
//                "/Volumes/igv/data/broad/hic/human/GSM455137_30305AAXX.1.maq.hic.summary.binned.txt.gz",
//                "/Volumes/igv/data/broad/hic/human/GSM455138_30305AAXX.2.maq.hic.summary.binned.txt.gz"};

        Map<String, Integer> chrIndexMap = new HashMap<String, Integer>();
        for (Chromosome chr : chromosomes) {
            if (chr != null && chr.getIndex() > 0)
                chrIndexMap.put(chr.getName(), chr.getIndex());
        }


        // Limit calcs to 10KB
        int[] gridSizeArray = new int[8];
        for (int i = 0; i < 8; i++) {
            gridSizeArray[i] = MainWindow.zoomBinSizes[i];
        }

        DensityCalculation[] calcs = new DensityCalculation[gridSizeArray.length];
        for (int z = 0; z < gridSizeArray.length; z++) {
            calcs[z] = new DensityCalculation(chromosomes, gridSizeArray[z]);
        }

        for (String path : paths) {
            AsciiPairIterator iter = new AsciiPairIterator(path);
            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                if (pair.getChr1().equals(pair.getChr2())) {
                    int dist = Math.abs(pair.getPos1() - pair.getPos2());


                    String chrName1 = pair.getChr1();
                    Integer index = chrIndexMap.get(chrName1);

                    if (index != null) {   // Make sure we know this chromosome
                        for (int z = 0; z < gridSizeArray.length; z++) {
                            calcs[z].addDistance(index, dist);
                        }
                    }
                }
            }
        }
        for (int z = 0; z < gridSizeArray.length; z++) {
            calcs[z].computeDensity();
        }

        outputDensities(calcs, new File("/xchip/igv/dev/hic/testFiles/test.hic.densities"));
    }


    private static void outputDensities(DensityCalculation[] calcs, File outputFile) throws IOException {

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
    }

    /**
     * Return a map of zoom level -> DensityFunction
     * @param is
     * @return
     * @throws IOException
     */
    public static Map<Integer, DensityFunction> readDensities(InputStream is) throws IOException {

        LittleEndianInputStream les = new LittleEndianInputStream(new BufferedInputStream(is));

        int nZooms = les.readInt();
        Map<Integer, DensityFunction> densityMap = new HashMap<Integer, DensityFunction>();
        // TODO -- Its assumed densities are in number order and indeces match resolutions.  This is fragile,
        // encode resolutions in the next round
        for (int i = 0; i < nZooms; i++) {

            DensityCalculation calc = new DensityCalculation(is);
            //   public DensityFunction(int gridSize, double[] densities, Map<Integer, Double> normFactors) {

            int gridSize = calc.getGridSize();
            double[] densities = calc.getDensityAvg();
            Map<Integer, Double> normFactors = calc.getNormalizationFactors();
            densityMap.put(i, new DensityFunction(gridSize, densities, normFactors));
        }

        return densityMap;

    }
}
