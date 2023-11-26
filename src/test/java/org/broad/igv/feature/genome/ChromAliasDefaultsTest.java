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
        ncibMockGenome.setChromAliasSource(new ChromAliasDefaults("GRCh*", ncbiNames));

        ucscNames = new ArrayList<>();
        for (int i = 1; i < 23; i++) {
            ucscNames.add("chr" + i);
        }
        ucscNames.add("chrX");
        ucscNames.add("chrY");
        ucscNames.add("chrM");
        AtomicInteger idx2 = new AtomicInteger();
        List<Chromosome> ucscChromosomes = ncbiNames.stream().map(name -> new Chromosome(idx2.getAndIncrement(), name, 0)).collect(Collectors.toList());
        ucscMockGenome = new Genome("hg*", ucscChromosomes);
        ucscMockGenome.setChromAliasSource(new ChromAliasDefaults("hg**", ucscNames));
    }

    @Test
    public void getCanonicalChromosomeName() {

        assertEquals(ncibMockGenome.getCanonicalChrName("chrX"), "23");
        assertEquals(ncibMockGenome.getCanonicalChrName("chrM"), "MT");

        assertEquals(ucscMockGenome.getCanonicalChrName("23"), "chrX");
        assertEquals(ucscMockGenome.getCanonicalChrName("MT"), "chrM");

    }


}