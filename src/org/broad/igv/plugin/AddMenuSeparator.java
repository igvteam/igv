/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin;

import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.ui.IGV;

import javax.swing.*;

/**
 * This is probably a case of over-design.
 * This plugin exists solely to structure the tools menu to our liking, it
 * just adds a separator to the menu.
 * @author jacob
 * @date 2013-Apr-12
 */
public class AddMenuSeparator implements IGVPlugin{
    @Override
    public void init() {
        IGV.getInstance().addOtherToolMenu(new JSeparator());
    }
}
