/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.scatterplot;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * User: jrobinso
 * Date: Jan 28, 2010
 */

public class ScatterData {

    private int count;
    private Map<String, double[]> dataMap;
    private Map<String, String[]> symbolMap;
    String[] samples;

    public ScatterData(int count, String[] samples, Map<String, double[]> dataMap, Map<String, String[]> symbolMap) {
        this.count = count;
        this.samples = samples;
        this.dataMap = dataMap;
        this.symbolMap = symbolMap;
    }

    public int getCount() {
        return count;
    }

    public List<String> getLabels() {
        return new ArrayList(dataMap.keySet());
    }

    public double[] getData(String label) {
        return dataMap.get(label);
    }

    public String[] getSymbols(String label) {
        return symbolMap.get(label);
    }

    public String[] getSamples() {
        return samples;
    }

    public static ScatterData getTestData() {

        int count = 300;
        double[] cnData = new double[count];
        double[] exprData = new double[count];
        double[] methData = new double[count];
        String[] treated = new String[count];
        String[] hyperMutated = new String[count];
        String[] samples = new String[count];

        for (int i = 0; i < count; i++) {
            samples[i] = "Sample " + i;
            cnData[i] = ((double) i / count) * 4 * (2 * Math.random());
            exprData[i] = -2 + ((double) i / count) * 4 * (2 * Math.random());
            methData[i] = 1 - i * Math.random() / count;
            treated[i] = Math.random() > 0.5 ? "AA" : "B";
            hyperMutated[i] = Math.random() > 0.9 ? "Y" : "N";

            System.out.println(cnData[i] + "\t" + exprData[i] + "\t" + methData[i]);
        }

        Map<String, double[]> dataMap = new LinkedHashMap();
        dataMap.put("Copy number", cnData);
        dataMap.put("Expression", exprData);
        dataMap.put("Methylation", methData);

        Map<String, String[]> symbolMap = new LinkedHashMap();
        symbolMap.put("Treated", treated);
        symbolMap.put("Hyper mutatated", hyperMutated);

        return new ScatterData(count, samples, dataMap, symbolMap);


    }

    public static void main(String[] args) {
        ScatterData sd = getTestData();

        for (String l : sd.getLabels()) {
            System.out.println(l);
        }
    }

}
