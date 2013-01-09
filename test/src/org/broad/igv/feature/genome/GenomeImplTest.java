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

package org.broad.igv.feature.genome;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class GenomeImplTest {
    @Test
    public void testGetNCBIName() throws Exception {

        String ncbiID = "gi|125745044|ref|NC_002229.3|";
        String ncbiName = "NC_002229.3";

        assertEquals(ncbiName, GenomeImpl.getNCBIName(ncbiID));

    }


    @Test
    public void testOrderChromos() throws Exception {
        //Scratch work for now, test methods for separating
        //contigs into "small" and "large"
        String indexPath = TestUtils.DATA_DIR + "fasta/CE.cns.all.fa.fai";
        Sequence seq = new MockSequence(indexPath);
        Genome genome = new GenomeImpl("GenomeImpleTest", "GenomeImplTest", seq);
        List<String> actNames = genome.getAllChromosomeNames();

        String[] expNames = {"chr1", "chr2", "chr3", "chrX", "C121713571", "scaffold22502"};
        int[] expInds = {0, 1, 2, 21, 22, actNames.size() - 1};
        int counter = 0;
        for (int expInd : expInds) {
            String expName = expNames[counter];
            String actName = actNames.get(expInd);
            assertEquals(expName, actName);
            counter++;
        }

    }

    /**
     * Class which loads FastaIndex and returns information contained therein,
     * but doesn't actually load full fasta file. For testing
     */
    private class MockSequence implements Sequence {

        private FastaIndex index;

        public MockSequence(String fastaIndexPath) throws IOException {
            this.index = new FastaIndex(fastaIndexPath);
        }

        @Override
        public byte[] getSequence(String chr, int start, int end) {
            return new byte[0];
        }

        @Override
        public byte getBase(String chr, int position) {
            return 0;
        }

        @Override
        public List<String> getChromosomeNames() {
            return new ArrayList<String>(index.getSequenceNames());
        }

        @Override
        public int getChromosomeLength(String chrname) {
            return index.getSequenceSize(chrname);
        }
    }
}
