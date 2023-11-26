package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class ChromAliasFileTest {

    private static Genome mockGenome;

    @org.junit.BeforeClass
    public static void setup() {
        List<Chromosome> chromosomeList = new ArrayList<>();
        chromosomeList.add(new Chromosome(0, "NC_007194.1", 0));
        mockGenome = new Genome("NC_007194.1", chromosomeList);
    }
    
    @Test
    public void getChromosomeName() throws IOException {
        String path = TestUtils.DATA_DIR +  "genomes/GCF_000002655.1.chromAlias.txt";
        ChromAliasFile chromAlias = new ChromAliasFile(path, mockGenome);
        assertEquals("NC_007194.1", chromAlias.getChromosomeName("CM000169.1")) ;
        assertEquals("NC_007194.1", chromAlias.getChromosomeName( "1"));
        assertEquals("NC_007194.1", chromAlias.getChromosomeName( "chr1"));
    }

    @Test
    public void search() throws IOException {

        String path = TestUtils.DATA_DIR +  "genomes/GCF_000002655.1.chromAlias.txt";
        ChromAliasFile chromAliasSource = new ChromAliasFile(path, mockGenome);
        ChromAlias chromAlias =  chromAliasSource.search("1");
        assertEquals(chromAlias.getChr(), "NC_007194.1");
        assertEquals(chromAlias.get("chr"), "NC_007194.1");
        assertEquals(chromAlias.get("genbank"), "CM000169.1") ;
        assertEquals(chromAlias.get("ncbi"), "1");
        assertEquals(chromAlias.get("ucsc"), "chr1");
    }
}
