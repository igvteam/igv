/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.exceptions;

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
