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

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


/**
 * @author jrobinso
 */
public class SegmentedBinaryLocalReader implements SegmentedBinaryReader {

    String filePath;

    List<String> chromosomes;

    List<String> sampleNames;

    Map<String, String> dataFileMap = new HashMap();

    Map<String, String> attrs;

    public SegmentedBinaryLocalReader(String filePath) {
        this.filePath = filePath;
    }

    /**
     * Return the set of sample names in this data set
     *
     * @return
     */
    public List<String> getSampleNames() {
        if (sampleNames == null) {
            sampleNames = readAsStrings("data/samples.txt");
        }
        return sampleNames;

    }

    public String getStringAttribute(String key) {
        return getAttributes().get(key);
    }

    public Map<String, String> getAttributes() {
        if (attrs == null) {
            attrs = new HashMap();
            List<String> lines = readAsStrings("data/attributes.txt");
            if (lines != null) {
                for (String kv : lines) {
                    String[] tokens = kv.split("=");
                    attrs.put(tokens[0], tokens[1]);
                }
            }
        }
        return attrs;

    }

    public SegmentedChromosomeData getChromosomeData(String chr) {

        return readChromosomeData(chr);

    }

    /**
     * Return the segment start locations for the sample and chromosome
     *
     * @param sample the sample name
     * @param chr    the chromosome name
     * @return array of end locations
     */
    public int[] getStartLocations(String sample, String chr) {
        String entryName = "data/" + sample + "/" + chr + "/starts.bin";
        return readAsInts(entryName);
    }

    /**
     * Return the segment end locations for the sample and chromosome
     *
     * @param sample the sample name
     * @param chr    the chromosome name
     * @return array of end locations
     */
    public int[] getEndLocations(String sample, String chr) {
        String entryName = "data/" + sample + "/" + chr + "/ends.bin";
        return readAsInts(entryName);
    }

    /**
     * Return the data values for the given sample and chromosom
     *
     * @param sample the sample name
     * @param chr    the chromosome name
     * @return data values
     */
    public float[] getValues(String sample, String chr) {
        String entryName = "data/" + sample + "/" + chr + "/values.bin";
        return readAsFloats(entryName);
    }

    private List<String> readAsStrings(String entryName) {
        List<String> strings = new ArrayList();
        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(new File(filePath));

            ZipEntry entry = zipFile.getEntry(entryName);
            if (entry == null) {
                return null;
            }

            BufferedReader br = new BufferedReader(new InputStreamReader(zipFile.getInputStream(entry)));
            String nextLine = null;
            while ((nextLine = br.readLine()) != null) {
                strings.add(nextLine.trim());
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                zipFile.close();
            } catch (IOException e) {
                // Close quietly
            }

        }
        return strings;

    }

    private SegmentedChromosomeData readChromosomeData(String chr) {
        ZipFile zipFile = null;
        BufferedInputStream is = null;
        try {
            String entryName = "data/" + chr + ".bin";
            zipFile = new ZipFile(new File(filePath));
            ZipEntry entry = zipFile.getEntry(entryName);
            if (entry == null) {
                return null;
            }

            is = new BufferedInputStream(zipFile.getInputStream(entry));

            SegmentedChromosomeData cd = new SegmentedChromosomeData();
            cd.deserialize(is);
            return cd;

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                zipFile.close();
            } catch (IOException e) {
                // Close quietly
            }
        }

    }

    /**
     * Read the contents of a zip entry as an array of ints.
     *
     * @param entryName name of the entry
     * @return array of ints
     */
    private int[] readAsInts(String entryName) {

        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(new File(filePath));
            ZipEntry entry = zipFile.getEntry(entryName);
            if (entry == null) {
                return null;
            }

            DataInputStream is = new DataInputStream(
                    new BufferedInputStream(zipFile.getInputStream(entry)));

            int nInts = (int) (entry.getSize() / 4);
            int[] ints = new int[nInts];
            for (int i = 0; i <
                    nInts; i++) {
                ints[i] = is.readInt();
            }

            return ints;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                zipFile.close();
            } catch (IOException e) {
                // close quietly
            }
        }

    }

    /**
     * Read the content of a zip entry as an array of floats.
     *
     * @param entryName name of the entry
     * @return contents as an array of floats
     */
    private float[] readAsFloats(String entryName) {

        ZipFile zipFile = null;
        try {
            zipFile = new ZipFile(new File(filePath));
            ZipEntry entry = zipFile.getEntry(entryName);
            if (entry == null) {
                return null;
            }

            DataInputStream is = new DataInputStream(
                    new BufferedInputStream(zipFile.getInputStream(entry)));

            int nInts = (int) (entry.getSize() / 4);
            float[] values = new float[nInts];
            for (int i = 0; i <
                    nInts; i++) {
                values[i] = is.readFloat();
            }

            return values;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            try {
                zipFile.close();
            } catch (IOException e) {

            }
        }


    }

}


