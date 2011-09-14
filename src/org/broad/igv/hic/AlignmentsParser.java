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

package org.broad.igv.hic;

import org.broad.igv.hic.data.Chromosome;
import org.broad.igv.hic.data.Dataset;
import org.broad.igv.hic.data.Matrix;
import org.broad.tribble.util.ParsingUtils;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * New format (Aug 2011)
 *         EAS1533_0010:4:77:10734:2666#0 X 8057 + SI-EAS1533_0010:4:77:10734:2666#0/2 X 5 -
 * @author jrobinso
 * @date Aug 2, 2010
 */
public class AlignmentsParser {

    public static Map<String, Integer> chromOrdinals = new HashMap<String, Integer>();

    static {
        chromOrdinals.put("2L", 0);
        chromOrdinals.put("2R", 1);
        chromOrdinals.put("3L", 2);
        chromOrdinals.put("3R", 3);
        chromOrdinals.put("4", 4);
        chromOrdinals.put("X", 5);
        chromOrdinals.put("U", 6);
    }

    public static Matrix readMatrix(InputStream is, String genomeId, int c1, int c2) throws IOException {
        String[] tokens = new String[9];
        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        Matrix matrix = new Matrix(c1, c2);
        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            int nTokens = ParsingUtils.split(nextLine, tokens, ' ');
            if (nTokens <= 9) {
                int chr1 = 0;
                int chr2 = 0;
                try {
                    // TODO -- use gnomeID to look up chromosome table
                    chr1 = chromOrdinals.get(tokens[1]);
                    chr2 = chromOrdinals.get(tokens[5]);
                } catch (Exception e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

                if ((c1 == chr1 && c2 == chr2) || (c1 == chr2 && c2 == chr1)) {
                    int pos1 = Integer.parseInt(tokens[2]);
                    int pos2 = Integer.parseInt(tokens[6]);
                    incrementCount(matrix, chr1, pos1, chr2, pos2);
                }
            }
        }
        matrix.parsingComplete();
        return matrix;
    }


    private static void incrementCount(Matrix matrix, int chr1, int pos1, int chr2, int pos2) {

        if (chr2 > chr1) {
            //transpose
            int tc2 = chr2;
            int tp2 = pos2;
            chr2 = chr1;
            pos2 = pos1;
            chr1 = tc2;
            pos1 = tp2;
        }

        matrix.incrementCount(pos1, pos2);
    }


}
