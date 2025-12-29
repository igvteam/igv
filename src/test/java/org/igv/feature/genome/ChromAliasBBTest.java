package org.igv.feature.genome;

import org.igv.feature.Chromosome;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class ChromAliasBBTest {
    private static Genome mockGenome;
    @org.junit.BeforeClass
    public static void setup() {
        List<Chromosome> chromosomeList = new ArrayList<>();
        chromosomeList.add(new Chromosome(0, "NC_007194.1", 0));
        mockGenome = new Genome("NC_007194.1", chromosomeList);
    }

    @Test
    public void getChromosomeName() throws IOException {
        String path = "test/data/genomes/GCF_000002655.1.chromAlias.bb";
        ChromAliasSource chromAlias = new ChromAliasBB(path, mockGenome);
        assertEquals("NC_007194.1", chromAlias.search("CM000169.1").getChr());
        assertEquals("NC_007194.1", chromAlias.search( "1").getChr());
        assertEquals("NC_007194.1", chromAlias.search( "chr1").getChr());
    }

    @Test
    public void search() throws IOException {

        String path = "test/data/genomes/GCF_000002655.1.chromAlias.bb";
        ChromAliasSource chromAliasSource = new ChromAliasBB(path, mockGenome);
        ChromAlias chromAlias =  chromAliasSource.search("1");
        assertEquals(chromAlias.get("genbank"), "CM000169.1") ;
        assertEquals(chromAlias.get("ncbi"), "1");
        assertEquals(chromAlias.get("ucsc"), "chr1");
    }

}