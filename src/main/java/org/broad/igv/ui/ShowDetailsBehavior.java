package org.broad.igv.ui;

/**
 * Created by jrobinso on 7/6/17.
 */
public enum ShowDetailsBehavior {

    HOVER("Show Details on Hover"),
    CLICK("Show Details on Click"),
    NEVER("Never Show Details");

    private final String label;

    private ShowDetailsBehavior(String label) {
        this.label = label;
    }

    public String getLabel() {
        return this.label;
    }
}
