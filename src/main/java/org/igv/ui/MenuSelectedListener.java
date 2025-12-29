package org.igv.ui;

import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

/**
 * A menu listener which provides default noop implementations of menuDeselected and menuCancelled
 */
public interface MenuSelectedListener extends MenuListener {

    @Override
    void menuSelected(MenuEvent e);

    @Override
    default void menuDeselected(MenuEvent e) {}

    @Override
    default void menuCanceled(MenuEvent e) {}
}
