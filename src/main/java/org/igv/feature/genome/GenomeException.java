package org.igv.feature.genome;


/**
 * @author eflakes
 */
public class GenomeException extends RuntimeException {

    public GenomeException(String message) {
        super(message);
    }

    public GenomeException(String message, Throwable e) {
        super(message, e);
    }
}