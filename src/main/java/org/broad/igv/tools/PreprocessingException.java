/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.tools;

/**
 * @author jrobinso
 */
public class PreprocessingException extends RuntimeException {

    public PreprocessingException(String msg) {
        super(msg);
    }

    PreprocessingException(String message, Throwable cause) {
        super(message, cause);
    }
}
