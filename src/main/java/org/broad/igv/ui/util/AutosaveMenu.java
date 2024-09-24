package org.broad.igv.ui.util;

import org.broad.igv.session.autosave.SessionAutosaveManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.MenuSelectedListener;
import org.broad.igv.ui.action.OpenSessionMenuAction;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Optional;

/**
 * A menu for listing currently available autosaved sessions that can be opened
 */
public class AutosaveMenu extends JMenu {

    public AutosaveMenu() {
        this("Autosaved Sessions");
    }

    private AutosaveMenu(String name) {
        super(name);

        this.addMenuListener((MenuSelectedListener) e -> fillAutosaveList());
    }

    /**
     * Fills this menu with options for opening each of the files in the autosave directory
     */
    private void fillAutosaveList() {
        // Get the exit session autosave file
        Optional<File> exitAutosave = SessionAutosaveManager.getExitSessionAutosaveFile();
        // Get a current list of timed autosaves
        File[] timedAutosaves = SessionAutosaveManager.getTimedSessionAutosaveFiles();
        Arrays.sort(timedAutosaves, Collections.reverseOrder());
        // Remove what's in the menu right now
        this.removeAll();
        // Add the exit session autosave file if it exists
        if(exitAutosave.isPresent()) {
            add(MenuAndToolbarUtils.createMenuItem(
                    new OpenSessionMenuAction(exitAutosave.get().getAbsolutePath(), IGV.getInstance())
            ));
            // Separate this from the timed autosaves
            if(timedAutosaves.length > 0) {
                addSeparator();
            }
        }
        // Create a menu item for each of the timed autosave files and add it to the menu
        for (File timedAutosave : timedAutosaves) {
            add(MenuAndToolbarUtils.createMenuItem(
                    new OpenSessionMenuAction(timedAutosave.getAbsolutePath(), IGV.getInstance())
            ));
        }
    }
}
