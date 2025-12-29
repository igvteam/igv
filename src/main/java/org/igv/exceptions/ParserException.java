package org.igv.exceptions;

/**
 * Created by IntelliJ IDEA.
 * User: nazaire
 * Date: Jul 13, 2009
 */
public class ParserException extends RuntimeException {
    private long lineNumber = -1;
    private String line;

    public ParserException(String message, long lineNumber) {
        super(message);

        setLineNumber(lineNumber);
    }

    public ParserException(String message, long lineNumber, String line) {
        super(message);

        setLineNumber(lineNumber);

        setLine(line);
    }

    public ParserException(String message, Throwable th, long lineNumber) {
        super(message, th);

        setLineNumber(lineNumber);
    }

    public ParserException(String message, Throwable th, long lineNumber, String line) {
        super(message, th);

        setLineNumber(lineNumber);

        setLine(line);
    }


    public void setLineNumber(long lineNumber) {
        this.lineNumber = lineNumber;
    }


    public void setLine(String line) {
        if (line != null) {
            this.line = line;
            if (line.length() > 500) {
                this.line = line.substring(0, 500);
            }
        }
    }

    public String getMessage() {
        String message = super.getMessage();
        if (message == null)
            message = "";

        if (getCause() != null) {
            if (line != null) {
                return "Failed to parse line " + lineNumber + ":\n"
                        + "Cause\n  " + getCause().getClass().getSimpleName() + ": " + message;

            } else {
                return "Failed to parse line " + lineNumber + ":\n"
                        + "\t" + line + "\n"
                        + "Cause\n  " + getCause().getClass().getSimpleName() + ": " + message;
            }
        }

        if (line != null) {
            return "Failed to parse line " + lineNumber + ":\n"
                    + "\t" + line + "\n"
                    + "\n  " + message;
        }

        return "Failed to parse line " + lineNumber + ":\n"
                + "\n  " + message;
    }
}
