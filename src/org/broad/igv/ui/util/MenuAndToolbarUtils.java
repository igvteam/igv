/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
