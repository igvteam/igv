
// should perhaps be in igv.data or some other location
package org.broad.igv.feature.basepair;

import java.awt.*; // does this have Color definition?
import java.util.List;
import java.util.HashMap;

public class BasePairData{
    // should have enumerated colors here
    public Color[] colors;

    // should hold list of helices with same color, and enum ID for arc color
    public int[][] startLeftNucs;
    public int[][] startRightNucs;
    public int[][] endLeftNucs;
    public int[][] endRightNucs;

    // convert lists to simple arrays
    public BasePairData(List<Color> inputColors, HashMap<Color,List<Object>> rowsByColor){
        colors = new Color[inputColors.size()];
        inputColors.toArray(colors);
        startLeftNucs = new int[colors.length][0];
        startRightNucs = new int[colors.length][0];
        endLeftNucs = new int[colors.length][0];
        endRightNucs = new int[colors.length][0];
        for (int i=0 ; i< colors.length ; i++){
            Color color = colors[i];
            System.out.println(color);
            List<Object> rows = rowsByColor.get(color);

            startLeftNucs[i] = new int[rows.size()];
            startRightNucs[i] = new int[rows.size()];
            endLeftNucs[i] = new int[rows.size()];
            endRightNucs[i] = new int[rows.size()];

            for (int j=0 ; j<rows.size() ; j++){
                // this is ugly
                List<Object> row = (List<Object>) rows.get(j);
                startLeftNucs[i][j] = (Integer) (row.get(0));
                startRightNucs[i][j] = (Integer) (row.get(1));
                endLeftNucs[i][j] = (Integer) (row.get(2));
                endRightNucs[i][j] = (Integer) (row.get(3));
                //System.out.format("%d\t%d\t%d\t%d\n",startLeftNucs[i][j],startRightNucs[i][j],endLeftNucs[i][j],endRightNucs[i][j]);
            }
        }
        //System.out.println("Colors: ");
        //for (int i=0; i<colors.length; i++){
        //    System.out.println(colors[i]);
        //}
    }
}