package org.broad.igv.feature.tribble;

/**
 * Thrown when an index is not found, but required.  This class extends Exception, as opposed to RuntimeException, so
 * clients will be forced to deal with it.
 *
 * @author jrobinso
 *         Date: 11/16/13
 *         Time: 10:10 PM
 */
public class TribbleIndexNotFoundException extends Exception{

    public TribbleIndexNotFoundException(String message) {
        super(message);
    }
}
