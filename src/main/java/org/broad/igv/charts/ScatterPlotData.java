/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.charts;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;


/*
*   Container class to encapsulate ScatterPlotFrame data for display purposes
*
*   Note: the data reference is copied for peformance reasons, but the data
*   is not protected from modification. This could be considerd a benefit
*   or a defect.
*
*   @author martind
* */
public class ScatterPlotData {

    String title;

    /**
     * Array of sample names
     */
    private String[] sampleNames;

    /**
     * Map of data type (symbol) =>  data array.  Data order must coincide with sample name array
     */
    private Map<String, double[]> dataMap;  // Data measurement keyname, data values

    /**
     * Sample attribute map.   Attribute heading => Array of values.  Value order must coincide with sample name array
     */
    private Map<String, String[]> symbolValueMap; // Attribute keyname, symbol map

    private int [] mutationCount;


    private int sampleCount;

    // Cached map key lists
    private ArrayList<String> dataTypes;   // list of names used for data map keys
    private ArrayList<String> symbolNames;   // list of names used for symbol map keys


    /*
   *    Constructor for containing IGV scatterplot data.
   *
   *    Parameters:
   *    sampleNames: data sample names for N samples
   *    symbolMap: symbol names for data sample attributes; e.g. "treated",
   *            "hyper mutated", etc.
   *            Each symbol map has an array of N symbol entries.
   *    dataMap: data values for data measurements; e.g. "copy number" ,
   *            "expression", "methylation", etc.
   *            Each data map has an array of N data entries.
   *
   * */
    public ScatterPlotData(String title,
                           String[] sampleNames,
                           Map<String, String[]> symbolMap,
                           Map<String, double[]> dataMap,
                           int [] mutationCount) {
        this.title = title;
        this.sampleNames = sampleNames;
        this.sampleCount = sampleNames.length;
        this.symbolValueMap = symbolMap;
        this.symbolNames = new ArrayList<String>(symbolValueMap.keySet());
        this.dataMap = dataMap;
        this.dataTypes = new ArrayList<String>(this.dataMap.keySet());
        this.mutationCount = mutationCount;

        // TODO -- validation

    }


    public String getTitle() {
        return title;
    }

    /*
   *   Returns the number of  data samples
   * */
    public int getSampleCount() {
        return sampleCount;
    }

    /*
    *   Returns array of data sample names
    * */
    public String[] getSampleNames() {
        return sampleNames;
    }

    public List<String> getCategories() {
        return symbolNames;
    }

    /*
    *   Returns the number of symbol keys in the SymbolMap
    * */
    public int getDataMapSize() {
        return dataMap.size();
    }

    /*
    *   Returns array of  data names
    * */
    public ArrayList<String> getDataNames() {
        return dataTypes;
    }

    /*
   *   Returns data measurement key name for given data map index
   * */
    public String getDataKeyName(int index) {
        return dataTypes.get(index);
    }

    /*
    *   Returns a map index for a given data measurement key name
    * */
    public int getDataKeyIndex(String keyName) {

        int labelIndex = -1;    // not assigned yet
        int nLabels = dataTypes.size();

        for (int index = 0; index < nLabels; ++index) {
            if (dataTypes.get(index).equals(keyName)) {
                labelIndex = index;
                // only finds the first instance
                break;
            }
        }
        return labelIndex;
    }

    /*
    *   Returns sample data measurement for a given data measurement key name
    *   and sample index.
    * */
    public double getDataKeyValue(String keyName, int sampleIndex) {
        return dataMap.get(keyName)[sampleIndex];
    }

    /*
   *   Returns an array of data measurements for a given data name key
   * */
    public double[] getDataValues(String name) {
        return dataMap.get(name);
    }

    /*
    *   Returns the number of symbol keys in the SymbolMap
    * */
    public int getSymbolMapSize() {
        return symbolValueMap.size();
    }

    /*
   *   Returns array of  symbol names
   * */
    public ArrayList<String> getSymbolNames() {
        return symbolNames;
    }

    /*
   *   Returns symbol key name for given map index
   * */
    public String getSymbolKeyName(int index) {
        return symbolNames.get(index);
    }

    /*
    *   Returns a map index for a given symbol key name
    * */
    public int getSymbolKeyIndex(String name) {

        int symbolIndex = -1;    // not assigned yet
        int nSymbols = symbolNames.size();

        for (int index = 0; index < nSymbols; ++index) {
            if (symbolNames.get(index).equals(name)) {
                symbolIndex = index;
                break;
            }
        }
        return symbolIndex;
    }

    /*
    *   Returns a sample values for a symbol key name and sample index.
    * */
    public String getSymbolKeyValue(String keyName, int sampleIndex) {
        return symbolValueMap.get(keyName)[sampleIndex];
    }

    /*
    *   Returns an array of all sample values for a symbol key name.
    * */
    public String[] getSymbolValues(String keyName) {
        return symbolValueMap.get(keyName);
    }

    /*
    *   Returns an array of unique named values found
    *   for a symbol attribute key.
    * */
    public String[] getAttributeCategories(String attribute) {

        // get the symbol values and load the first unique variant
        String[] symbolValues = getSymbolValues(attribute);
        ArrayList<String> uniqueVals = new ArrayList<String>();

        // the "no-value" value
        uniqueVals.add("");

        // find all unique symbol values for attribute key
        for (String x : symbolValues) {
            if (x != null && !uniqueVals.contains(x)) {
                uniqueVals.add(x);
            }
        }

        int nCategories = uniqueVals.size();
        String[] categories = new String[nCategories];
        uniqueVals.toArray(categories);

        return categories;
    }

    /*
     * @deprecated
    *   Returns a sample description String for a particular
    *   data sample built from all sample point information.
    **/
    public String getSampleDescription(int sampleIndex, boolean isHTML) {
        String description = "<html>" + sampleNames[sampleIndex] + ": <br>";
        String keyName;
        int keyCount;
        int index;

        Double Value;
        keyCount = getDataMapSize();
        for (index = 0; index < keyCount; ++index) {
            keyName = getDataKeyName(index);
            Value = getDataKeyValue(keyName, sampleIndex);
            description += keyName + " = " + Value.toString() + "<br>";
        }

        String Sample;
        keyCount = getSymbolMapSize();

        if (isHTML) {
            for (index = 0; index < keyCount; ++index) {
                keyName = getSymbolKeyName(index);
                Sample = getSymbolKeyValue(keyName, sampleIndex);
                description += keyName + " = " + Sample + "<br>";
            }

            // end of description
            description += "</html>";
        } else {
            for (index = 0; index < keyCount; ++index) {
                keyName = getSymbolKeyName(index);
                Sample = getSymbolKeyValue(keyName, sampleIndex);
                if (index < keyCount - 1)
                    description += keyName + " = " + Sample + ", ";
                else
                    description += keyName + " = " + Sample;
            }

        }
        return description;
    }

    /**
     * Marker for mutated samples
     */
    public int[] getMutationCount() {
        return mutationCount;
    }
}
