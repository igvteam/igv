package org.igv.ui.dnd;

import java.awt.event.MouseAdapter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GhostDropAdapter extends MouseAdapter {
    protected GhostGlassPane glassPane;
    protected String action;

    private List listeners;

    public GhostDropAdapter(GhostGlassPane glassPane, String action) {
        this.glassPane = glassPane;
        this.action = action;
        this.listeners = new ArrayList();
    }

    public void addGhostDropListener(GhostDropListener listener) {
        if (listener != null)
            listeners.add(listener);
    }

    public void removeGhostDropListener(GhostDropListener listener) {
        if (listener != null)
            listeners.remove(listener);
    }

    protected void fireGhostDropEvent(GhostDropEvent evt) {
        Iterator it = listeners.iterator();
        while (it.hasNext()) {
            ((GhostDropListener) it.next()).ghostDropped(evt);
        }
    }
}
