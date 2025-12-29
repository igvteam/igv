package org.broad.igv.sam.reader;

/**
 * @author jrobinso
 */
public class UnsortedFileException extends RuntimeException {

    public UnsortedFileException(Throwable cause) {
        super(cause);
    }

    public UnsortedFileException(String message, Throwable cause) {
        super(message, cause);
    }

    public UnsortedFileException(String message) {
        super(message);
    }

    public UnsortedFileException() {
    }
}

