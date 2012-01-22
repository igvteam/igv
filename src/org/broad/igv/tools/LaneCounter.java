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

package org.broad.igv.tools;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.reader.AlignmentQueryReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Feb 1, 2011
 */


public class LaneCounter {


    //int intervalStart = 152668596;
    //int intervalEnd = 154746948;
    //static String alignmentFile = "http://iwww.broadinstitute.org/igvdata/sigma/aggregation/G8556/564/v2/564.bam";


    /**
     * arg[0] => genome id
     * arg[2] => comma delimited list of bam files
     * arg[3] =>  query interval string
     */
    public static void main(String[] args) throws IOException {

        Globals.setHeadless(true);

        Genome genome = IgvTools.loadGenome(args[0], false);
        if (genome == null) {
            throw new PreprocessingException("Genome could not be loaded: hg18");
        }

        run(genome, args[1], args[2]);

    }

    public static void run(Genome genome, String alignmentFileList, String queryInterval) {

        String[] tokens = queryInterval.split(":");
        String chr = tokens[0];
        String[] tokens2 = tokens[1].split("-");
        int intervalStart = Integer.parseInt(tokens2[0]);
        int intervalEnd = Integer.parseInt(tokens2[1]);

        byte[] ref = genome.getSequence(chr, intervalStart, intervalEnd);
        if(ref == null) {
            return;
        }


        String[] tmp = alignmentFileList.split(",");
        for (String alignmentFile : tmp) {

            System.out.println(alignmentFile + "   " + queryInterval);

            AlignmentQueryReader reader = null;
            CloseableIterator<Alignment> iter = null;

            try {

                reader = AlignmentReaderFactory.getReader(alignmentFile, true);
                iter = reader.query(chr, intervalStart, intervalEnd, false);
                Map<String, ReadGroupCount> counts = new HashMap();

                while (iter != null && iter.hasNext()) {
                    Alignment alignment = iter.next();
                    if (passFilter(alignment)) {

                        if (alignment.getMappingQuality() == 0) {
                            // TODO -- mq zero event
                        }

                        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                        if (blocks != null) {
                            int lastBlockEnd = -1;
                            for (AlignmentBlock block : blocks) {

                                if (!block.isSoftClipped()) {

                                    byte[] bases = block.getBases();
                                    int blockStart = block.getStart();
                                    int adjustedStart = block.getStart();
                                    int adjustedEnd = block.getEnd();


                                    for (int pos = adjustedStart; pos < adjustedEnd; pos++) {

                                        if (pos < intervalStart) {
                                            continue;
                                        }
                                        if (pos >= intervalEnd) {
                                            break;
                                        }
                                        int baseIdx = pos - blockStart;
                                        if (bases != null && baseIdx >= 0 && baseIdx < bases.length) {
                                            String rg = alignment.getReadGroup();

                                            ReadGroupCount rgc = counts.get(rg);
                                            if (rgc == null) {
                                                rgc = new ReadGroupCount(rg);
                                                counts.put(rg, rgc);
                                            }
                                            rgc.totalBases++;

                                            byte base = bases[baseIdx];
                                            byte refBase = ref[pos - intervalStart];
                                            if (base != refBase) {
                                                rgc.mismatches++;
                                            }
                                        }
                                    }

                                    lastBlockEnd = block.getEnd();
                                }
                            }
                        }
                    }
                }

                outputReadCounts(counts);
            }

            catch (Exception e) {
                e.printStackTrace();
            }

            finally {
                iter.close();
                try {
                    reader.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
    }

    private static void outputReadCounts(Map<String, ReadGroupCount> counts) {
        for (ReadGroupCount rgc : counts.values()) {
            float percent = (100.0f * rgc.mismatches) / rgc.totalBases;
            System.out.println(rgc.rg + "\t" + rgc.totalBases + "\t" + rgc.mismatches + "\t" + percent);
        }
    }

    private static boolean passFilter(Alignment alignment) {
        return alignment.isMapped() &&
                !alignment.isDuplicate() &&
                alignment.getMappingQuality() > 0 &&
                !alignment.isVendorFailedRead();
    }


    static class ReadGroupCount {

        String rg;
        int totalBases = 0;
        int mismatches = 0;

        ReadGroupCount(String rg) {
            this.rg = rg;
        }

    }


}
