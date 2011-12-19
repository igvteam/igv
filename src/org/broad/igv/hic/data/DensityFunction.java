package org.broad.igv.hic.data;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 12/5/11
 */
public class DensityFunction {


    //String densityFile = "HumanAugust.densities.txt";
    int[] positions;
    double[] density;
    int gridSize;
    private int nPoints;
    Map<Integer, Double> normFactors;

    public DensityFunction(int gridSize, double[] densities, Map<Integer, Double> normFactors) {
        this.gridSize = gridSize;
        this.density = densities;
        this.nPoints = densities.length;
        this.normFactors = normFactors;
    }

    public double getDensity(int chrIdx, int distance, int binSize) {

        double ratio = ((double) binSize) / gridSize;
        double mult = 1.0; // ratio*ratio;


        double normFactor = normFactors.containsKey(chrIdx) ? normFactors.get(chrIdx) : 1.0;
        mult *= normFactor;

        int grid = distance / gridSize;
        if (grid >= nPoints) {
            return density[nPoints - 1] * mult;
        } else {
            return density[grid] * mult;
        }
    }


//    private void init() {
//
//        try {
//            InputStream is = DensityFunction.class.getResourceAsStream(densityFile);
//            BufferedReader br = new BufferedReader(new InputStreamReader(is));
//            nPoints = Integer.parseInt(br.readLine());
//            gridSize = Integer.parseInt(br.readLine());
//
//            int idx = 0;
//            positions = new int[nPoints];
//            density = new double[nPoints];
//            String nextLine;
//            while ((nextLine = br.readLine()) != null && idx < nPoints) {
//                String[] tokens = nextLine.split("\t");
//                positions[idx] = Integer.parseInt(tokens[0]);
//                density[idx] = Double.parseDouble(tokens[1]);
//                idx++;
//
//            }
//        } catch (IOException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }
//    }
//

}
