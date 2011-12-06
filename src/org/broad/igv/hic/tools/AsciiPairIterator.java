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

package org.broad.igv.hic.tools;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public class AsciiPairIterator implements PairIterator {

    AlignmentPair nextPair = null;
    BufferedReader reader;


    public AsciiPairIterator(String path) throws IOException {
        this.reader = org.broad.igv.util.ParsingUtils.openBufferedReader(path);
        advance();
    }

    private void advance() {

        try {
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {

                String[] tokens = nextLine.split("\\s+");
                int nTokens = tokens.length;
                if (nTokens <= 9) {
                    String chrom1 = tokens[1];
                    String chrom2 = tokens[5];
                    int pos1 = Integer.parseInt(tokens[2]);
                    int pos2 = Integer.parseInt(tokens[6]);

                    nextPair = new AlignmentPair(chrom1, pos1, chrom2, pos2);
                    return;
                }

            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        nextPair = null;

    }

    public boolean hasNext() {
        return nextPair != null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public AlignmentPair next() {
        AlignmentPair p = nextPair;
        advance();
        return p;
    }

    public void remove() {
        // Not implemented
    }

    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

}
