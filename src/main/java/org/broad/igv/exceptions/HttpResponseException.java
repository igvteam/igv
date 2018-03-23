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

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Class for converting HTTP error response codes to exceptions
 *
 * @author Jim Robinson
 * @date Jul 27, 2011
 */
public class HttpResponseException extends IOException {

    int statusCode;
    String message;
    String details;

    public HttpResponseException(int statusCode, String message, String details) {
        this.statusCode = statusCode;
        this.details = details;
        if(message != null && message.length() > 0) {
            this.message = message;
        }
        else {
            switch (statusCode) {
                case 407:
                    this.message = "Proxy authentication required (status code " + statusCode + ")";
                case 403:
                    this.message = "Access Forbidden (status code " + statusCode + ")";
                case 404:
                    this.message = "File not found (status code " + statusCode + ")";
                case 401:
                    this.message = "Not authorized (status code " + statusCode + ")";
                default:
                    this.message = "HTTP access error (status code " + statusCode + ")";
            }
        }
    }

    public int getStatusCode() {
        return statusCode;
    }

    @Override
    public String getMessage() {

        if(this.details == null) {
            return message;
        }
        else {
            return message + " <br> " + details;
        }


    }
}
