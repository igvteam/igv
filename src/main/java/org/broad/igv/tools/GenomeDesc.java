/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tools;

import java.util.List;

/**
 * @author jrobinso
 */
public class GenomeDesc {

    String genomeId;
    List<Chrom> chromosomes;
    long length;

    public long getLength() {
        return length;
    }


    //public long getCumulativeOffset(chr)
    //genome.getChromosomeLength(chr)
    class Chrom {

        String name;
        int length;
    }
}
