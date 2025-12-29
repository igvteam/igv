package org.broad.igv.sam.reader;

/**
 * @author jrobinso
 * @date Oct 20, 2010
 */
public class IndexNotFoundException extends RuntimeException {

    String samFile;
    public IndexNotFoundException(String samFile) {
        this.samFile = samFile;
    }

    public String getSamFile() {
        return samFile;
    }

    public String getMessage() {
        return "Could not find index for file: " + samFile;
    }
}
