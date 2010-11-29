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

package org.broad.igv.feature;

import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

/**
 * 1. matches - Number of bases that match that aren't repeats
 * 2. misMatches - Number of bases that don't match
 * 3. repMatches - Number of bases that match but are part of repeats
 * 4. nCount - Number of 'N' bases
 * 5. qNumInsert - Number of inserts in query
 * 6. qBaseInsert - Number of bases inserted in query
 * 7. tNumInsert - Number of inserts in target
 * 8. tBaseInsert - Number of bases inserted in target
 * 9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
 * 10. qName - Query sequence name
 * 11. qSize - Query sequence size
 * 12. qStart - Alignment start position in query
 * 13. qEnd - Alignment end position in query
 * 14. tName - Target sequence name
 * 15. tSize - Target sequence size
 * 16. tStart - Alignment start position in target
 * 17. tEnd - Alignment end position in target
 * 18. blockCount - Number of blocks in the alignment (a block contains no gaps)
 * 19. blockSizes - Comma-separated list of sizes of each block
 * 20. qStarts - Comma-separated list of starting positions of each block in query
 * 21. tStarts - Comma-separated list of starting positions of each block in target
 * <p/>
 * Example:
 * track name=fishBlats description="Fish BLAT" useScore=1
 * 59 9 0 0 1 823 1 96 +- FS_CONTIG_48080_1 1955 171 1062 chr22 47748585 13073589 13073753 2 48,20,  171,1042,  34674832,34674976,
 * 59 7 0 0 1 55 1 55 +- FS_CONTIG_26780_1 2825 2456 2577 chr22 47748585 13073626 13073747 2 21,45,  2456,2532,  34674838,34674914,
 * 59 7 0 0 1 55 1 55 -+ FS_CONTIG_26780_1 2825 2455 2576 chr22 47748585 13073627 13073748 2 45,21,  249,349,  13073627,13073727,
 *
 * @author jrobinso
 * @date Aug 5, 2010
 */
public class PSLParser extends UCSCParser {

    @Override
    public boolean isFeatureFile(ResourceLocator locator) {
        String path = locator.getPath();
        return path.endsWith(".psl") || path.endsWith(".psl.gz")  ||
                path.endsWith(".pslx") || path.endsWith(".pslx.gz");
    }

    @Override
    protected BasicFeature parseLine(String[] tokens, int nTokens) {

        if (nTokens < 21) {
            // log.info("Skipping line ")
            return null;
        }
        int tSize = Integer.parseInt(tokens[14]);
        String chr = tokens[13];
        int start = Integer.parseInt(tokens[15]); // IS PSL 1 or ZERO based,  closed or open?

        String strandString = tokens[8];
        Strand strand = strandString.startsWith("+") ? Strand.POSITIVE : Strand.NEGATIVE;

        boolean gNeg = false;
        if (strandString.length() > 1) {
            gNeg = strandString.charAt(1) == '-';
        }

        BasicFeature f = new BasicFeature();
        f.setName(tokens[9]);
        f.setChr(chr);
        f.setStart(start);
        f.setEnd(Integer.parseInt(tokens[16]));
        f.setStrand(strand);

        int exonCount = Integer.parseInt(tokens[17]);
        String[] exonSizes = tokens[18].split(",");
        String[] startsBuffer = tokens[20].split(",");

        if (startsBuffer.length == exonSizes.length && exonCount == startsBuffer.length) {
            for (int i = 0; i < startsBuffer.length; i++) {
                int exonSize = Integer.parseInt(exonSizes[i]);
                int exonStart = Integer.parseInt(startsBuffer[i]);
                if (gNeg) {
                    exonStart = tSize - exonStart - exonSize;
                }
                int exonEnd = exonStart + exonSize;
                Exon exon = new Exon(chr, exonStart, exonEnd, strand);
                f.addExon(exon);
            }
        } else {
            // TODO -- warn
        }


        return f;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
