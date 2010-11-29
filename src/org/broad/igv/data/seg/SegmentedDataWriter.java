/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;

import java.io.*;
import java.util.*;
import java.util.zip.ZipOutputStream;

/**
 * @author jrobinso
 */
public class SegmentedDataWriter {

    private static Logger log = Logger.getLogger(SegmentedDataWriter.class);

    private boolean compress = true;

    SegmentedAsciiDataSet dataset;

    String rootNodeName;

    File outputFile;

    TrackType trackType;

    public SegmentedDataWriter(SegmentedAsciiDataSet dataset, File outputFile) {
        this(dataset, outputFile, false, TrackType.COPY_NUMBER);
    }

    public SegmentedDataWriter(SegmentedAsciiDataSet dataset, File outputFile, boolean compress, TrackType type) {
        this.dataset = dataset;
        this.compress = compress;
        this.outputFile = outputFile;
        this.rootNodeName = "data";
        this.trackType = type;
    }

    public void writeContents() {

        ZipOutputStream zos = null;

        try {
            zos = new ZipOutputStream(new FileOutputStream(outputFile));
            List<String> samples = dataset.getDataHeadings();
            String[] sampleArray = samples.toArray(new String[]{});
            Set<String> chromosomes = dataset.getChromosomes();

            List<String> attributes = new ArrayList();
            attributes.add("type=" + trackType.toString());
            attributes.add("logNormalized=" + dataset.isLogNormalized());
            writeStrings("attributes.txt", attributes, zos);

            writeStrings("samples.txt", samples, zos);
            writeStrings("chromosomes.txt", chromosomes, zos);


            ArrayList<String> tmp = new ArrayList(dataset.getChromosomes());
            tmp.add(Globals.CHR_ALL);
            for (String chr : tmp) {
                Map<String, int[]> starts = new HashMap();
                Map<String, int[]> ends = new HashMap();
                Map<String, float[]> values = new HashMap();
                for (String sample : samples) {

                    List<LocusScore> segments = chr.equals(Globals.CHR_ALL) ? 
                            dataset.getWholeGenomeScores(sample) :
                            dataset.getSegments(sample, chr);
                    int nPts = segments == null ? 0 : segments.size();
                    int[] s = new int[nPts];
                    int[] e = new int[nPts];
                    float[] v = new float[nPts];
                    for (int i = 0; i < nPts; i++) {
                        s[i] = segments.get(i).getStart();
                        e[i] = segments.get(i).getEnd();
                        v[i] = segments.get(i).getScore();
                    }

                    starts.put(sample, s);
                    ends.put(sample, e);
                    values.put(sample, v);

                    //writeSegments(sample, chr, segments, zos);
                }
                SegmentedChromosomeData cd = new SegmentedChromosomeData(
                        sampleArray, starts, ends, values);
                writeChromsomeData(chr, cd, zos);

            }

        } catch (IOException e) {
            e.printStackTrace();

        } finally {
            try {
                zos.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private void writeStrings(String filename, Collection<String> strings, ZipOutputStream zos) {
        try {
            // Put strings in a byte buffer.  We need to know the total size in bytes
            StringBuffer buff = new StringBuffer();
            for (String s : strings) {
                buff.append(s);
                buff.append('\n');
            }
            byte[] bytes = buff.toString().getBytes();

            String entryName = rootNodeName + "/" + filename;
            Utilities.createZipEntry(entryName, bytes, zos, compress);

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
        }
    }

    /**
     * Write the segment data from a single chromosome for a sample.  The data
     * is organized in the following hierarchical structure:
     * <p/>
     * data/<sample>/<chr>/starts.bin
     * data/<sample>/<chr>/ends.bin
     * data/<sample>/<chr>/values.bin
     *

     */
    private void writeChromsomeData(String chr, SegmentedChromosomeData cd, ZipOutputStream zos) {
        try {

            ByteArrayOutputStream bytes = new ByteArrayOutputStream(1000);
            BufferedOutputStream bos = new BufferedOutputStream(bytes);
            cd.serialize(bytes);
            bos.close();

            byte[] byteArray = bytes.toByteArray();

            String entryName = rootNodeName + "/" + chr + ".bin";
            Utilities.createZipEntry(entryName, byteArray, zos, compress);

        } catch (IOException ex) {
            log.error("Error writing segments", ex);
        } finally {
        }
    }

    public static void main(String[] args) {


        String inputFile = null;
        String outputFile = null;
        String genomeId = null;
        boolean compress = true;
        TrackType trackType = TrackType.COPY_NUMBER;
        if (args.length > 2) {
            inputFile = args[0].trim();
            outputFile = args[1].trim();
            genomeId = args[2].trim();
        } else {
            System.out.println("Arguments: inputFile  outputFile  genomeId  [track type (optional)");
            System.exit(-1);
        }
        if (args.length > 3) {
            try {
                trackType = TrackType.valueOf(args[3].trim());
            } catch (Exception e) {
                System.out.println("Unknown track type: " + args[3]);
            }
        }

        System.out.println("Track type=" + trackType.toString());
        //inputFile = "/Users/jrobinso/IGVTestData/Demo/TCGA/broad.seg";
        //outputFile = "test.seg.zip";

        // TODO -- remove the need for this hack
        System.out.println("Setting genome: " + genomeId);
       GenomeManager.getInstance().setGenomeId(genomeId);

        SegmentedAsciiDataSet ds = new SegmentedAsciiDataSet(new ResourceLocator(inputFile));

        SegmentedDataWriter writer = new SegmentedDataWriter(ds, new File(outputFile), compress, trackType);

        writer.writeContents();

    }

    public void setCompress(

            boolean compress) {
        this.compress = compress;
    }
}
