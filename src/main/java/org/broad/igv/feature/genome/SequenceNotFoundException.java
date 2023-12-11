package org.broad.igv.feature.genome;

public class SequenceNotFoundException extends RuntimeException {
    public SequenceNotFoundException(String message) {
        super(message);
    }
}
