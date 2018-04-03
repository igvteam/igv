package org.broad.igv.event;

/**
 * Created by jrobinso on 2/26/17.
 *
 * Represents view shift left or right,  e.g. moving via drag or arrow key.
 */
public class ShiftEvent {

    public int delta;

    public ShiftEvent(int delta) {
        this.delta = delta;
    }
}
