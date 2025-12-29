/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.tdf;

/**
 * @author jrobinso
 */
public interface IBFBin {

    public int getSize();

    public int getStartPosition(int idx);

    public int getEndPosition(int idx);

    public float getValue(int row, int idx);

}
