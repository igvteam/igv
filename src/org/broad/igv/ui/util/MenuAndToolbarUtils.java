/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * MenuAndToolbarUtils.java
 *
 * Created on November 7, 2007, 1:34 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.ui.util;

import org.broad.igv.ui.action.MenuAction;

import javax.swing.*;
import java.util.List;

/**
 * @author eflakes
 */
public class MenuAndToolbarUtils {


    static public JMenu createMenu(List<JComponent> menuItems, MenuAction action) {

        final JMenu menu = new JMenu();
        menu.setAction(action);

        for (JComponent menuItem : menuItems) {
            menu.add(menuItem);
        }

        String text = action.getToolTipText();
        if (text != null) {
            menu.setToolTipText(text);
        }

        return menu;
    }

    static public JMenuItem createMenuItem(MenuAction menuItemAction) {

        JMenuItem menuItem = new JMenuItem();
        menuItem.setAction(menuItemAction);

        String text = menuItemAction.getToolTipText();
        if (text != null) {
            menuItem.setToolTipText(text);
        }

        return menuItem;
    }

    static public JCheckBoxMenuItem createMenuItem(MenuAction menuItemAction, boolean isSelected) {

        JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
        menuItem.setSelected(isSelected);
        menuItem.setAction(menuItemAction);

        String text = menuItemAction.getToolTipText();
        if (text != null) {
            menuItem.setToolTipText(text);
        }

        return menuItem;
    }
}
