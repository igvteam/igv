/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui;

import org.broad.igv.dev.api.api;

/**
 * Identifier for each track panel which is displayed
 * in the IGV window
 * User: jacob
 * Date: 2013-Jan-31
 */
@api
public enum PanelName {
    FEATURE_PANEL(IGV.FEATURE_PANEL_NAME),
    DATA_PANEL(IGV.DATA_PANEL_NAME);

    private final String panelName;

    private PanelName(String panelName){
        this.panelName = panelName;
    }

    public String getName(){
        return this.panelName;
    }


}
