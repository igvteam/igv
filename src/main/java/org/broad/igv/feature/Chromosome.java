/*
Chromosome.java
 *
Created on June 29, 2007, 8:53 AM
 *
To change this template, choose Tools | Template Manager
and open the template in the editor.
 */
package org.broad.igv.feature;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * Simple representation of a chromosome.  Basically a name, length, and optionally a list of cytobands.
 *
 * @author jrobinso
 */
public class Chromosome {
    private String name;
    private  transient int index;  // Order in the chromosome (for convenience)
    private int length = 0;

    public Chromosome(int index, String name, int length) {
        this.index = index;
        this.name = name;
        this.length = length;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int ii) {
        this.index = ii;
    }

    /**
     * /**
     * Return the length of the chromosome
     */
    public int getLength() {
        return length;
    }

    public String getName() {
        return name;
    }

    public String toString() {
        return name;
    }

    @Override
    public boolean equals(Object obj) {
        return (obj instanceof Chromosome) &&
                ((Chromosome)obj).name.equals(name)
                && ((Chromosome)obj).length == length;
    }

    @Override
    public int hashCode() {
        return Objects.hash(name, length);
    }
}
