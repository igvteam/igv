
package org.igv.track;

import java.io.FileNotFoundException;

/**
 * @author jrobinso
 */
public interface TrackReader {

    TrackSet getTracks() throws FileNotFoundException;

}
