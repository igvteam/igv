
package org.igv.tools;

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
