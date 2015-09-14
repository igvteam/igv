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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ArrayHeapIntSorter;
import org.broad.igv.util.IntComparator;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * @author jrobinso
 */
public class WiggleDataset implements Dataset {

    Genome genome;
    private String name;
    private TrackProperties trackProperties;
    Map<String, IntArrayList> startLocationsMap = new HashMap();
    Map<String, IntArrayList> endLocationsMap = new HashMap();
    Map<String, FloatArrayList> dataMap = new HashMap();

    float dataMin = 0;
    float dataMax = 0;
    float percent90 = 0;
    float percent10 = 0;

    private Map<String, Integer> longestFeatureMap;
    private TrackType type = TrackType.OTHER;

    public WiggleDataset(Genome genome, String name) {
        this.genome = genome;
        this.name = name;
        this.trackProperties = new TrackProperties();
    }

    public void sort(Set<String> unsortedChromosomes) {
        for (String c : unsortedChromosomes) {
            String chr = genome.getChromosomeAlias(c);

            final IntArrayList starts = startLocationsMap.get(chr);
            int sz = starts.size();

            int[] indices = new int[sz];
            for (int i = 0; i < indices.length; i++) {
                indices[i] = i;
            }

            (new ArrayHeapIntSorter()).sort(indices, new IntComparator() {

                public int compare(int arg0, int arg1) {
                    return starts.get(arg0) - starts.get(arg1);
                }
            });


            int[] sortedStarts = reorder(indices, startLocationsMap.get(chr));
            int[] sortedEnds = reorder(indices, endLocationsMap.get(chr));
            float[] sortedData = reorder(indices, dataMap.get(chr));

            startLocationsMap.put(chr, new IntArrayList(sortedStarts));
            endLocationsMap.put(chr, new IntArrayList(sortedEnds));
            dataMap.put(chr, new FloatArrayList(sortedData));
        }

    }

    private float[] reorder(int[] indices, FloatArrayList values) {
        int size = values.size();
        if (indices.length != size) {
            throw new IllegalArgumentException(
                    "Index array length not equal to size");
        }
        float[] reorderedValues = new float[size];
        for (int i = 0; i < size; i++) {
            reorderedValues[i] = values.get(indices[i]);
        }
        return reorderedValues;
    }

    private int[] reorder(int[] indices, IntArrayList values) {
        int size = values.size();
        if (indices.length != size) {
            throw new IllegalArgumentException(
                    "Index array length not equal to size");
        }
        int[] reorderedValues = new int[size];
        for (int i = 0; i < size; i++) {
            reorderedValues[i] = values.get(indices[i]);
        }
        return reorderedValues;
    }


    public void addDataChunk(String chr, IntArrayList starts, IntArrayList ends, FloatArrayList data) {
        IntArrayList startLocations = this.startLocationsMap.get(chr);
        if (startLocations == null) {
            this.startLocationsMap.put(chr, starts);
        } else {
            startLocations.addAll(starts);
        }

        if (ends != null) {
            IntArrayList endLocations = this.endLocationsMap.get(chr);
            if (endLocations == null) {
                this.endLocationsMap.put(chr, ends);
            } else {
                endLocations.addAll(ends);
            }
        }

        FloatArrayList dataArray = this.dataMap.get(chr);
        if (dataArray == null) {
            this.dataMap.put(chr, data);
        } else {

            dataArray.addAll(data);
        }
        float[] d = data.toArray();
        for (int i = 0; i < d.length; i++) {
            dataMax = Math.max(dataMax, d[i]);
            dataMin = Math.min(dataMin, d[i]);
        }


    }


    public float getDataMin() {
        return dataMin;
    }


    public float getDataMax() {
        return dataMax;
    }


    public String getName() {
        return name;
    }


    public TrackType getType() {
        return type;
    }

    public String[] getChromosomes() {
        return startLocationsMap.keySet().toArray(new String[]{});
    }

    boolean containsChromosome(String chr){
        String datChr = this.genome != null ? genome.getChromosomeAlias(chr) : chr;
        return startLocationsMap.containsKey(datChr);
    }


    public String[] getTrackNames() {
        return new String[]{getName()};
    }


    public int[] getStartLocations(String chr) {
        IntArrayList startLocations = this.startLocationsMap.get(chr);
        if (startLocations == null) {
            return null;
        } else {
            return startLocations.toArray();
        }
    }

    public int[] getEndLocations(String chr) {
        IntArrayList endLocations = this.endLocationsMap.get(chr);
        if (endLocations == null) {
            return null;
        } else {
            return endLocations.toArray();
        }
    }

    public float[] getData(String heading, String chr) {
        FloatArrayList data = this.dataMap.get(chr);
        if (data == null) {
            return null;
        } else {
            return data.toArray();
        }
    }

    public String[] getFeatureNames(String chr) {
        return null;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public void setName(String name) {
        this.name = name;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public Integer getLongestFeature(String chr) {
        return longestFeatureMap == null ? 1000 :
                longestFeatureMap.containsKey(chr) ? longestFeatureMap.get(chr) : 1;
    }

    public void setLongestFeatureMap(Map<String, Integer> longestFeatureMap) {
        this.longestFeatureMap = longestFeatureMap;
    }

    public void setType(TrackType type) {
        this.type = type;
    }

    public void setPercent10(float percent10) {
        this.percent10 = percent10;
    }

    public float getPercent10() {
        return percent10;
    }

    public void setPercent90(float percent90) {
        this.percent90 = percent90;
    }

    public float getPercent90() {
        return percent90;
    }
}
