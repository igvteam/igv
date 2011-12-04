package org.broad.igv.ui.event;

import java.util.EventObject;

/**
 * @author Jim Robinson
 * @date 12/2/11
 */
public class AlignmentTrackEvent extends EventObject {

    public enum Type {SPLICE_JUNCTION, VISIBILITY_WINDOW, RELOAD}

    private Type type;

    public AlignmentTrackEvent(Object source, Type type) {
        super(source);
        this.type = type;
    }

    public Type getType() {
        return type;
    }
}
