package org.igv.feature;

import org.igv.feature.genome.GenomeException;


/**
 * @author eflakes
 */
public class MaximumContigGenomeException extends GenomeException {

    public MaximumContigGenomeException(String message) {
        super(message);
    }

    public MaximumContigGenomeException(String message, Throwable e) {
        super(message, e);
    }
}