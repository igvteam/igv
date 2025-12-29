package org.igv.exceptions;

/**
 * Author: nazaire
 * Date: Jul 8, 2009
 */
public class DataLoadException extends RuntimeException {

    private String message;
    private String fileName;

    public DataLoadException(String message) {
        this.message = message.replace("<html>", "");
    }

    public DataLoadException(String message, String fileName) {
        if(message != null) this.message = message.replace("<html>", "");
        this.fileName = fileName;
    }

    public String getMessage() {
        return fileName == null ? message : "An error occurred while accessing:    " + fileName + "<br>" + message;
    }
}