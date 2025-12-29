/*
 * MenuAndToolbarUtils.java
 *
 * Created on November 7, 2007, 1:34 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.igv.ui.util;

import org.igv.ui.action.MenuAction;

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
