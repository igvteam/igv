/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
