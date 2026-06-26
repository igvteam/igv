package org.igv.circview.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A genome assembly: an ordered list of chromosomes drawn around the circle.
 *
 * <p>Mirrors the {name, id, chromosomes} input to setAssembly() in
 * circularView.js. Chromosome names are shortened (leading "chr" stripped) and
 * a color is assigned from {@link org.igv.circview.util.ChrColors} when absent.
 */
public final class Assembly {

    private final String name;
    private final String id;
    private final List<Chromosome> chromosomes;

    public Assembly(String name, String id, List<Chromosome> chromosomes) {
        this.name = name;
        this.id = id;
        this.chromosomes = new ArrayList<>(chromosomes);
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public List<Chromosome> getChromosomes() {
        return chromosomes;
    }

    /** Total length in base pairs across all chromosomes. */
    public long totalBp() {
        long total = 0;
        for (Chromosome c : chromosomes) {
            total += c.getBpLength();
        }
        return total;
    }

    public Chromosome getChromosome(String name) {
        for (Chromosome c : chromosomes) {
            if (c.getName().equals(name)) {
                return c;
            }
        }
        return null;
    }
}
