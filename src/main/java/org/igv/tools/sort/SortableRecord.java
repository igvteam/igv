package org.igv.tools.sort;

/**
 * @author mnazaire
 */
public class SortableRecord {
    private String chromosome;
    private int start;
    private String text;

    public SortableRecord(String chromosome, int start, String text) {
        this.chromosome = chromosome;
        this.start = start;
        this.text = text;
    }

    public String getChromosome() {
        return chromosome;
    }

    public int getStart() {
        return start;
    }

    public String getText() {
        return text;
    }
}
