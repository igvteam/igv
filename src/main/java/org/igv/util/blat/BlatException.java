package org.igv.util.blat;

public class BlatException extends RuntimeException {

    public BlatException(String message) {
        super(message);
    }

    public BlatException(String message, Throwable cause) {
        super(message, cause);
    }
}
