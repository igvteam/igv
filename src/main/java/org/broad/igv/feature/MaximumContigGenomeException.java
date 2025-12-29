package org.broad.igv.feature;

import org.broad.igv.feature.genome.GenomeException;


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