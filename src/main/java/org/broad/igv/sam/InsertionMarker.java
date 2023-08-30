package org.broad.igv.sam;

import java.util.Objects;

/**
 * Created by jrobinso on 1/13/17.
 */
public class InsertionMarker {

    public int position;
    public int size;

    public InsertionMarker(int position, int size) {
        this.position = position;
        this.size = size;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        InsertionMarker that = (InsertionMarker) o;
        return position == that.position && size == that.size;
    }

    @Override
    public int hashCode() {
        return Objects.hash(position, size);
    }
}
