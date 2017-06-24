package org.broad.igv.feature.genome.fasta;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by jrobinso on 6/23/17.
 */
public class FastaBlockCompressedSequenceTest {


    @Test
    public void findBlockContaining() throws Exception {

        String fasta = "/Users/jrobinso/Downloads/hg38/hg38.fa.gz";

        FastaBlockCompressedSequence seq = new FastaBlockCompressedSequence(fasta);

        FastaBlockCompressedSequence.Mapping mapping = seq.findBlockContaining(1000000);

        assertNotNull(mapping);

    }

}
