package org.broad.igv.feature.genome.fasta;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by jrobinso on 6/23/17.
 */
public class FastaBlockCompressedSequenceTest {


    @Test
    public void findBlockContaining() throws Exception {

        String fasta = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.gz";

        FastaBlockCompressedSequence seq = new FastaBlockCompressedSequence(fasta);

        FastaBlockCompressedSequence.Mapping mapping = seq.findBlockContaining(1000000);

        assertNotNull(mapping);

    }

}
