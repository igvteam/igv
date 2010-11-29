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

package org.broad.igv.goby;

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.SamQueryReaderFactory;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jul 12, 2010
 * Time: 11:25:52 AM
 * To change this template use File | Settings | File Templates.
 */
public class GobyAlignmentQueryReaderTest {

    String thmFile = "test/data/goby/DLTTEJH-Bullard-HBR-SRR037439.tmh";
    GobyAlignmentQueryReader reader;

    @Before
    public void setup() throws IOException {
        reader = new GobyAlignmentQueryReader(thmFile);
    }

    @Test
    public void testGetSequenceNames() throws Exception {

        Set<String> expectedSequences = new HashSet(Arrays.asList("NT_113900", "22", "NT_113910", "NT_113911", "MT",
                "3", "NT_113904", "2", "NT_113903", "1", "NT_113902", "NT_113901", "NT_113908", "7", "NT_113907", "6",
                "5", "NT_113906", "4", "NT_113905", "9", "NT_113909", "8", "19", "c5_H2", "17", "NT_113920", "M",
                "18", "NT_113921", "15", "16", "13", "14", "11", "12", "NT_113917", "21", "NT_113916", "20",
                "NT_113919", "NT_113918", "NT_113913", "NT_113912", "NT_113915", "NT_113914", "Y", "X", "NT_113932",
                "NT_113933", "NT_113930", "NT_113931", "NT_113929", "NT_113928", "10", "NT_113927", "NT_113926",
                "NT_113925", "NT_113924", "NT_113923", "NT_113898", "NT_113899", "NT_113949", "NT_113945", "NT_113890",
                "NT_113946", "NT_113947", "NT_113948", "NT_113953", "NT_113952", "NT_113955", "c6_QBL", "NT_113954",
                "NT_113951", "NT_113950", "NT_113884", "NT_113885", "NT_113882", "NT_113883", "NT_113888", "NT_113889",
                "NT_113886", "NT_113887", "NT_113938", "NT_113939", "NT_113936", "NT_113880", "NT_113937", "NT_113881",
                "NT_113934", "NT_113935", "NT_113944", "NT_113943", "NT_113942", "NT_113941", "NT_113940", "c22_H2",
                "NT_113879", "NT_113875", "NT_113876", "NT_113877", "NT_113878", "NT_113871", "NT_113872", "NT_113873",
                "NT_113874", "NT_113870", "NT_113958", "c6_COX", "NT_113956", "NT_113957", "NT_113962", "NT_113961",
                "NT_113960", "NT_113966", "NT_113965", "NT_113964", "NT_113963"));

        Set<String> seqs = reader.getSequenceNames();
        assertEquals(expectedSequences.size(), seqs.size());
        for (String s : seqs) {
            assertTrue(expectedSequences.contains(s));
        }
    }

    @Test
    public void testIterator() throws Exception {

        CloseableIterator<Alignment> iter =  reader.iterator();
        while(iter.hasNext()) {
            Alignment al = iter.next();
            System.out.println(al.getChr() + ":" + al.getStart() + "-" + al.getEnd());
        }
    }

    @Test
    public void testQuery() throws Exception {
    }
 
}
