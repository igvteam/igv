package org.igv.circview.model;

import java.awt.Color;

/**
 * A single reference chromosome (or contig) in an {@link Assembly}.
 *
 * <p>Mirrors the {name, bpLength, color} chromosome objects in the JS assembly
 * input (see hg19.js / setAssembly in circularView.js).
 */
public final class Chromosome {

    private final String name;
    private final long bpLength;
    private final Color color;

    public Chromosome(String name, long bpLength, Color color) {
        this.name = name;
        this.bpLength = bpLength;
        this.color = color;
    }

    public String getName() {
        return name;
    }

    public long getBpLength() {
        return bpLength;
    }

    public Color getColor() {
        return color;
    }

    /** Strip a leading "chr" prefix, matching shortChrName() in circularView.js. */
    public static String shortChrName(String chrName) {
        return chrName.startsWith("chr") ? chrName.substring(3) : chrName;
    }
}
