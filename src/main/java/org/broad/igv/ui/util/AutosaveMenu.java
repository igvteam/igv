package org.broad.igv.ui.util;

import org.broad.igv.session.autosave.SessionAutosaveManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.OpenSessionMenuAction;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import java.io.File;
import java.util.Arrays;

/**
 * A menu for listing currently available autosaved sessions that can be opened
 */
public class AutosaveMenu extends JMenu {

    public AutosaveMenu() {
        this("Autosaved Sessions");
    }

    private AutosaveMenu(String name) {
        super(name);

        this.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                fillAutosaveList();
            }

            @Override
            public void menuDeselected(MenuEvent e) {

            }

            @Override
            public void menuCanceled(MenuEvent e) {

            }
        });
    }

    private void fillAutosaveList() {
        // Get a current list of autosaves
        File[] autosaves = SessionAutosaveManager.getSessionAutosaveFiles();
        Arrays.sort(autosaves);
        // Remove what's in the menu right now
        this.removeAll();
        // Create a menu item for each and add it to the menu
        for(int i = 0; i < autosaves.length; i++) {
            add(MenuAndToolbarUtils.createMenuItem(
                    new OpenSessionMenuAction(autosaves[i].getAbsolutePath(), IGV.getInstance())
            ));
        }
    }
}
