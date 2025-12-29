package org.igv.exceptions;

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
