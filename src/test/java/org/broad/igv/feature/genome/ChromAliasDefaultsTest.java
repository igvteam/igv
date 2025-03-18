package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class ChromAliasDefaultsTest {

    private static List<String> ucscNames;
    private static List<String> ncbiNames;

    private static Genome ucscMockGenome;

    private static Genome ncibMockGenome;

    @org.junit.BeforeClass
    public static void setup() {
        ncbiNames = new ArrayList<>();
        for (int i = 1; i < 25; i++) {
            ncbiNames.add(String.valueOf(i));
        }
        ncbiNames.add("MT");

        AtomicInteger idx = new AtomicInteger();
        List<Chromosome> ncbiChromosomes = ncbiNames.stream().map(name -> new Chromosome(idx.getAndIncrement(), name, 0)).collect(Collectors.toList());
        ncibMockGenome = new Genome("GRCh*", ncbiChromosomes);

        ucscNames = new ArrayList<>();
        for (int i = 1; i < 23; i++) {
            ucscNames.add("chr" + i);
        }
        ucscNames.add("chrX");
        ucscNames.add("chrY");
        ucscNames.add("chrM");
        AtomicInteger idx2 = new AtomicInteger();
        List<Chromosome> ucscChromosomes = ucscNames.stream().map(name -> new Chromosome(idx2.getAndIncrement(), name, 0)).collect(Collectors.toList());
        ucscMockGenome = new Genome("hg*", ucscChromosomes);
    }

    @Test
    public void getCanonicalChromosomeName() {

        assertEquals("23", ncibMockGenome.getCanonicalChrName("chrX"));
        assertEquals("23", ncibMockGenome.getCanonicalChrName("X"));
        assertEquals("MT", ncibMockGenome.getCanonicalChrName("chrM"));

        assertEquals("chrX", ucscMockGenome.getCanonicalChrName("23"));
        assertEquals("chrX", ucscMockGenome.getCanonicalChrName("X"));
        assertEquals("chrM", ucscMockGenome.getCanonicalChrName("MT"));

        assertEquals("chr1", ucscMockGenome.getCanonicalChrName("Chr1"));

    }


}