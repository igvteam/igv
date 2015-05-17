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
