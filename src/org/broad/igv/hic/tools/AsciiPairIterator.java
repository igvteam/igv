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


import org.broad.igv.Globals;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * @author Jim Robinson
 * @date 9/24/11
 */
public class AsciiPairIterator implements PairIterator {

    static Pattern WHITESPACE_PATTERN = Pattern.compile("\\s+");
    // Map of name -> index
    private Map<String, Integer> chromosomeOrdinals;
    AlignmentPair nextPair = null;
    BufferedReader reader;

    /**
     * A map of chromosome name -> chromosome string.  A private "intern" pool.  The java "intern" pool stores string
     * in perm space, which is rather limited and can cause us to run out of memory.
     */
    Map<String, String> stringInternPool = new HashMap();

    public AsciiPairIterator(String path, Map<String, Integer> chromosomeOrdinals) throws IOException {
        this.reader = org.broad.igv.util.ParsingUtils.openBufferedReader(path);
        this.chromosomeOrdinals = chromosomeOrdinals;
        advance();
    }

    /**
     * Read the next record
     * <p/>
     * D0J8AACXX120130:6:1101:1003:8700/1 15 61559113 0 D0J8AACXX120130:6:1101:1003:8700/2 15 61559309 16
     * D0J8AACXX120130:6:1101:1004:2368/1 10 26641879 16 D0J8AACXX120130:6:1101:1004:2368/2 9 12797549 0
     */
    private void advance() {

        try {
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
                int nTokens = tokens.length;
                if (nTokens == 6) {
                    String chrom1 = getInternedString(tokens[1]);
                    String chrom2 = getInternedString(tokens[4]);
                    if (chromosomeOrdinals.containsKey(chrom1) && chromosomeOrdinals.containsKey(chrom2)) {
                        int chr1 = chromosomeOrdinals.get(chrom1);
                        int chr2 = chromosomeOrdinals.get(chrom2);
                        int pos1 = Integer.parseInt(tokens[2]);
                        int pos2 = Integer.parseInt(tokens[5]);
                        nextPair = new AlignmentPair(chr1, pos1, chr2, pos2);
                        return;
                    }

                }
                else if (nTokens < 10) {
                    String chrom1 = getInternedString(tokens[1]);
                    String chrom2 = getInternedString(tokens[5]);

                    if (chromosomeOrdinals.containsKey(chrom1) && chromosomeOrdinals.containsKey(chrom2)) {
                        int chr1 = chromosomeOrdinals.get(chrom1);
                        int chr2 = chromosomeOrdinals.get(chrom2);
                        int pos1 = Integer.parseInt(tokens[2]);
                        int pos2 = Integer.parseInt(tokens[6]);
                        nextPair = new AlignmentPair(chr1, pos1, chr2, pos2);
                        return;
                    }
                }

            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        nextPair = null;

    }

    /**
     * Replace "aString" with a stored equivalent object, if it exists.  If it does not store it.  The purpose
     * of this class is to avoid running out of memory storing zillions of equivalent string.
     *
     * @param aString
     * @return
     */
    private String getInternedString(String aString) {
        String s = stringInternPool.get(aString);
        if (s == null) {
            s = new String(aString); // THe "new" will break any dependency on larger strings if this is a "substring"
            stringInternPool.put(aString, s);
        }
        return s;
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
