package org.igv.event;

import org.igv.ui.panel.ReferenceFrame;

/**
 * Events corresponding to a change in viewed area (chromosome, position, and/or zoom).
 * <p>
 * {@code Cause} derived events
 * should cause the data model (e.g. ReferenceFrame) to change their position, this will
 * typically be dispatched by a UI Component.
 * <p>
 * {@code Result} derived events should be dispatched after objects have changed
 * their position, typically to tell UI elements they should repaint
 * User: jacob
 * Date: 2013-Jan-30
 */
public final class ViewChange implements IGVEvent{

    public enum Type {ChromosomeChange, LocusChange}

    final boolean recordHistory;
    final public Type type;
    final public String chrName;
    final public double start;
    final public double end;
    final public ReferenceFrame referenceFrame;
    public boolean panning = false;
    public boolean fromPanning = false;

    private ViewChange(Type type, ReferenceFrame referenceFrame, String chrName, double start, double end, boolean recordHistory) {
        this.type = type;
        this.referenceFrame = referenceFrame;
        this.chrName = chrName;
        this.start = start;
        this.end = end;
        this.recordHistory = recordHistory;
    }


    public boolean recordHistory() {
        return this.recordHistory;
    }

    public static ViewChange ChromosomeChangeResult(ReferenceFrame referenceFrame, String chrName, boolean recordHistory) {
        return new ViewChange(Type.ChromosomeChange, referenceFrame, chrName, 0.0, 0.0, recordHistory);
    }

    public static ViewChange LocusChangeResult(ReferenceFrame referenceFrame, String chrName, double start, double end, boolean recordHistory) {
        return new ViewChange(Type.LocusChange, referenceFrame, chrName, start, end, recordHistory);
    }

    public static ViewChange LocusChangeResultPanning(ReferenceFrame referenceFrame, String chrName, double start, double end, boolean recordHistory) {
        ViewChange vc =  new ViewChange(Type.LocusChange, referenceFrame, chrName, start, end, recordHistory);
        vc.panning = true;
        return vc;
    }
}
