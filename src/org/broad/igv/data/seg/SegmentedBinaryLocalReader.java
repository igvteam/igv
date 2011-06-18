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


