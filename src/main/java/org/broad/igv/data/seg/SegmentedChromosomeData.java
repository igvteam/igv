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
package org.broad.igv.data.seg;

import org.apache.log4j.Logger;
import org.broad.igv.data.CharArrayList;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


/**
 * This class encapuslates data for a single chromosome from a segmented
 * dataset.   Its batch purpose is to provide a single object that can
 * easily be serialized and passed from server to IGV clients.
 *
 * @author jrobinso
 */
public class SegmentedChromosomeData {

    private static Logger log = Logger.getLogger(SegmentedChromosomeData.class);

    private String[] sampleNames;

    Map<String, int[]> startLocations;

    Map<String, int[]> endLocations;

    Map<String, float[]> valueMap;

    public SegmentedChromosomeData() {
    }

    public SegmentedChromosomeData(
            String[] sampleNames,
            Map<String, int[]> startLocations,
            Map<String, int[]> endLocations,
            Map<String, float[]> valueMap) {
        this.sampleNames = sampleNames;
        this.startLocations = startLocations;
        this.endLocations = endLocations;
        this.valueMap = valueMap;
    }

    /**
     * Serialize this instance to a stream.
     * <p/>
     * Format:
     * Sample names as a string :
     * sample1 tab sample2 tabnetc...  line-feed
     * Followed by (for each sample)
     * numberOfPoints (int)
     * start locations (int array)
     * end locations (int array)
     * values (float array)
     *
     * @param stream stream to serialize to
     */
    public void serialize(OutputStream stream) {

        try {
            DataOutputStream os = new DataOutputStream(new BufferedOutputStream(stream));

            for (String sample : getSampleNames()) {
                os.writeChars(sample);
                os.writeChar('\t');
            }
            os.writeChar('\n');
            for (String sample : getSampleNames()) {
                int[] starts = startLocations.get(sample);
                int[] ends = endLocations.get(sample);
                float[] values = valueMap.get(sample);
                int nPts = starts.length;

                assert ends.length == nPts;
                assert values.length == nPts;

                os.writeInt(nPts);
                for (int i = 0; i < nPts; i++) {
                    os.writeInt(starts[i]);
                }
                for (int i = 0; i < nPts; i++) {
                    os.writeInt(ends[i]);
                }
                for (int i = 0; i < nPts; i++) {
                    os.writeFloat(values[i]);
                }
            }
            os.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Initialize the object by deserializing the stream.
     *
     * @param stream
     * @see #serialize
     */
    public void deserialize(InputStream stream) {
        startLocations = new HashMap();
        endLocations = new HashMap();
        valueMap = new HashMap();
        try {

            DataInputStream is = new DataInputStream(new BufferedInputStream(stream));

            ArrayList<String> tmp = new ArrayList(100);
            char c = 0;
            try {
                c = is.readChar();
            } catch (EOFException e) {
                // This is expected when there is no data available for this chromosome.
                return;
            }
            while (c != '\n') {
                CharArrayList chars = new CharArrayList(10);
                while (c != '\t') {
                    chars.add(c);
                    c = is.readChar();
                }
                String nm = new String(chars.toArray()).trim();
                tmp.add(nm);
                c = is.readChar();
            }

            sampleNames = tmp.toArray(new String[]{});

            for (String sample : sampleNames) {
                int nPts = is.readInt();
                int[] starts = new int[nPts];
                for (int i = 0; i < nPts; i++) {
                    starts[i] = is.readInt();
                }
                int[] ends = new int[nPts];
                for (int i = 0; i < nPts; i++) {
                    ends[i] = is.readInt();
                }
                float[] values = new float[nPts];
                for (int i = 0; i < nPts; i++) {
                    values[i] = is.readFloat();
                }
                startLocations.put(sample, starts);
                endLocations.put(sample, ends);
                valueMap.put(sample, values);
            }

        }
        catch (IOException e) {
            log.error("Error loading segmented data", e);
        }
    }

    public String[] getSampleNames() {
        return sampleNames;
    }

    public int[] getStartLocations(String sampleName) {
        return startLocations.get(sampleName);
    }

    public int[] getEndLocations(String sampleName) {
        return endLocations.get(sampleName);
    }

    public float[] getValues(String sampleName) {
        return valueMap.get(sampleName);
    }
}
