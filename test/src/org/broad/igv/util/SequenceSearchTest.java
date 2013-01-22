/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util;

import org.broad.tribble.Feature;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2013-Jan-22
 */
public class SequenceSearchTest {


    String shortSeq = "GACTCTGACTGGACTCTGATCAG";

    @Test
    public void testBasicSearch() throws Exception{
        String motif = "TCTG";
        List<Feature> matches = SequenceSearch.search(motif, shortSeq.getBytes());

        assertEquals(2, matches.size());
        assertEquals(3, matches.get(0).getStart());
        assertEquals(7, matches.get(0).getEnd());

        assertEquals(14, matches.get(1).getStart());
        assertEquals(18, matches.get(1).getEnd());
    }

    @Test
    public void testConvertMotifToRegex_Basic() throws Exception{
        String motif = "ACTGACTGACTG";
        String regex = SequenceSearch.convertMotifToRegex(motif);
        assertEquals(motif, regex);
    }

    @Test
    public void testConvertMotifToRegex_02() throws Exception{
        String motif = "ACTGMACTGNACTSG";
        String regex = SequenceSearch.convertMotifToRegex(motif);

        assertTrue(regex.length() >= motif.length());
        assertEquals("ACTG[M,A,C]ACTG.ACT[S,G,C]G", regex);
    }
}
