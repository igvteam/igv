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
