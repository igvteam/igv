package org.igv.exceptions;

/**
 * User: nazaire
 */
public class ProbeMappingException extends RuntimeException {
    private String genome;
    private String fileName;

    public ProbeMappingException(String fileName, String genome) {
        super();
        this.fileName = fileName;
        this.genome = genome;
    }

    public String getGenome() {
        return genome;
    }

    public String getFileName() {
        return fileName;
    }

    public String getMessage() {
        return "The probes in the dataset could not be mapped to locations on the loaded genome (" + genome + ").";
    }
}
