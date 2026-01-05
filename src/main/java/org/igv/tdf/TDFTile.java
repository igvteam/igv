
package org.igv.tdf;

import java.io.IOException;

/**
 * @author jrobinso
 */
public interface TDFTile {
    int[] getStart();

    int[] getEnd();

    float[] getData(int trackNumber);

    String[] getNames();

    enum Type {
        fixedStep, variableStep, bed, bedWithName
    }

    ;

    public int getTileStart();

    public int getSize();

    public int getStartPosition(int idx);

    public int getEndPosition(int idx);

    public String getName(int idx);

    public float getValue(int row, int idx);

    public void writeTo(BufferedByteWriter fos) throws IOException;

}
