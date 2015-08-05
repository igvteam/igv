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

package org.broad.igv.tools.converters;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import org.apache.log4j.Logger;
import org.broad.igv.data.expression.GeneToLocusHelper;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.sort.SortableRecord;
import org.broad.igv.tools.sort.SortableRecordCodec;
import org.broad.igv.track.TrackType;

import java.io.*;
import java.util.Comparator;
import java.util.List;


/**
 * @author jrobinso
 * @date Oct 9, 2010
 * @deprecated Doesn't appear to be used anywhere, as of 23 May 2013
 */
@Deprecated
public class MageTabToIGVConverter {

    private static Logger log = Logger.getLogger(MageTabToIGVConverter.class);


    /**
     * Parse the file and output in ".igv" format
     *
     * @return
     */
    public static void convert(File inputFile, File outputFile, String probeResource,
                               int maxRecords, File tmpDir, Genome genome) throws IOException {

        GeneToLocusHelper locusHelper = new GeneToLocusHelper(probeResource);


        BufferedReader reader = null;
        PrintWriter writer = null;

        SortingCollection cltn = getSortingCollection(maxRecords, tmpDir);
        try {
            reader = new BufferedReader(new FileReader(inputFile));
            writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            // Parse the header line
            String headerLine = reader.readLine();
            String[] tokens = headerLine.split("\t");

            //The sample names in a GCT file start at column 2,
            int sampleStart = 2;


            String nextLine = null;
            TrackType dataType = TrackType.GENE_EXPRESSION;
            while ((nextLine = reader.readLine()) != null) {

                // A gct row can map to multiple loci, normally this indicates a problem with the probe
                DataRow row = new DataRow(nextLine);
                String probe = row.getProbe();
                if (probe.startsWith("cg")) {
                    dataType = TrackType.DNA_METHYLATION;
                }

                List<Locus> loci = locusHelper.getLoci(probe, row.getDescription(), genome.getId());
                if (loci == null || loci.isEmpty()) {
                    System.out.println("No locus found for: " + probe + "  " + row.getDescription());
                } else {
                    for (Locus locus : loci) {
                        String igvLine = locus.getChr() + "\t" + locus.getStart() + "\t" + locus.getEnd() + "\t" + probe +
                                row.getData();
                        cltn.add(new SortableRecord(locus.getChr(), locus.getStart(), igvLine));
                    }
                }
            }

            writer.println("#type=" + dataType.toString());
            writer.print("Chr\tStart\tEnd\tProbe");
            for (int i = sampleStart; i < tokens.length; i++) {
                writer.print("\t" + tokens[i]);
            }
            writer.println();
            // Ouputput the sorted file
            CloseableIterator<SortableRecord> iter = cltn.iterator();
            while (iter.hasNext()) {
                SortableRecord al = iter.next();
                writer.println(al.getText());

            }
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (writer != null) {
                writer.close();
            }
        }

    }


    static SortingCollection getSortingCollection(int maxRecords, File tmpDir) {


        SortableRecordCodec codec = new SortableRecordCodec();

        Comparator<SortableRecord> comp = new Comparator<SortableRecord>() {

            public int compare(SortableRecord o1, SortableRecord o2) {


                String chr1 = o1.getChromosome().replaceFirst("chr", "");
                String chr2 = o2.getChromosome().replaceFirst("chr", "");
                int s1 = Integer.MAX_VALUE;
                try {
                    s1 = Integer.parseInt(chr1);
                } catch (Exception e) {
                    // ignore
                }
                int s2 = Integer.MAX_VALUE;
                try {
                    s2 = Integer.parseInt(chr2);
                } catch (Exception e) {
                    // ignre
                }


                int t1 = s1 - s2;
                if (t1 == 0) {
                    chr1 = chr1.replace("M", "Z");
                    chr2 = chr2.replace("M", "Z");
                    t1 = chr1.compareTo(chr2);
                }
                if (t1 == 0) {
                    return (int) (o1.getStart() - o2.getStart());
                } else {
                    return t1;
                }
            }
        };


        return SortingCollection.newInstance(SortableRecord.class, codec, comp, maxRecords, tmpDir);
    }


    /**
     * Represents a row of data from a GCT or mage-tab file.  Using this class if more effecient than tokeninzing
     * the entire line.  Some GCT files have over a thousand columns and we're only interested in the first 2
     */
    static class DataRow {

        private String probe;
        private String description;
        private String data;

        DataRow(String string) {

            int firstTab = string.indexOf('\t');
            int secondTab = string.indexOf('\t', firstTab + 1);

            // TODO -- if either of the indeces firstTab or secondTab are < 0 throw an exception
            probe = string.substring(0, firstTab);
            description = string.substring(firstTab, secondTab);
            data = string.substring(secondTab);


        }

        private String getProbe() {
            return probe;
        }

        public String getDescription() {
            return description;
        }

        public String getData() {
            return data;
        }
    }

}
