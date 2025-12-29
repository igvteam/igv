/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.track;

import java.io.FileNotFoundException;

/**
 * @author jrobinso
 */
public interface TrackReader {

    TrackSet getTracks() throws FileNotFoundException;

}
