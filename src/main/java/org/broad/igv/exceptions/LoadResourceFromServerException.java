package org.broad.igv.exceptions;

/**
 * User: nazaire
 */
public class LoadResourceFromServerException extends RuntimeException {
    String resourceLocation = null;
    String exceptionClass = null;

    public LoadResourceFromServerException(String message) {
        super(message);
    }

    public LoadResourceFromServerException(String message, String resourceLocation) {
        super(message);
        this.resourceLocation = resourceLocation;
    }

    public LoadResourceFromServerException(String message, String resourceLocation, String exceptionClass) {
        super(message);
        this.resourceLocation = resourceLocation;
        this.exceptionClass = exceptionClass;
    }

    public String getMessage() {
        if (resourceLocation == null) {
            return "Unable to load a resource " +
                    "from the server \nCause\n " + super.getMessage();
        } else if (exceptionClass != null) {
            return "Unable to load the resource \n" +
                    "  " + resourceLocation +
                    "\nCause\n   " + exceptionClass + ": " + super.getMessage();
        } else {
            return "Unable to load the resource \n" +
                    resourceLocation +
                    " : \n" + super.getMessage();
        }
    }
}
