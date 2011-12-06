package org.broad.igv.hic.data;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * @author Jim Robinson
 * @date 12/5/11
 */
public class DensityFunction {


    String densityFile = "HumanAugust.densities.txt";
    int[] positions;
    double[] density;
    int gridSize;
    private int nPoints;

    public DensityFunction() {
        init();

    }

    public double getDensity(int distance) {
        int grid = distance / gridSize;
        if(grid >= nPoints) {
            return density[nPoints-1];
        }
        else {
            return density[grid];
        }
    }


    private void init() {

        try {
            InputStream is = DensityFunction.class.getResourceAsStream(densityFile);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            nPoints = Integer.parseInt(br.readLine());
            gridSize = Integer.parseInt(br.readLine());

            int idx = 0;
            positions = new int[nPoints];
            density = new double[nPoints];
            String nextLine;
            while ((nextLine = br.readLine()) != null && idx < nPoints) {
                String[] tokens = nextLine.split("\t");
                positions[idx] = Integer.parseInt(tokens[0]);
                density[idx] = Double.parseDouble(tokens[1]);
                idx++;

            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


}
