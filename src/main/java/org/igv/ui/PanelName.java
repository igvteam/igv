package org.igv.ui;

/**
 * Identifier for each track panel which is displayed
 * in the IGV window
 * @author jacob
 * @date 2013-Jan-31
 * @api
 */
public enum PanelName {
    ANNOTATION_PANEL(IGV.FEATURE_PANEL_NAME),
    DATA_PANEL(IGV.DATA_PANEL_NAME);

    private final String panelName;

    private PanelName(String panelName){
        this.panelName = panelName;
    }

    public String getName(){
        return this.panelName;
    }


}
