package org.broad.igv.sam;

/**
 * Created by jrobinso on 1/13/17.
 */
public class InsertionMarker {
    public int position;
    public int size;
    public int pixelPosition = -1;

    public InsertionMarker(int position, int size) {
        this.position = position;
        this.size = size;
    }
}
