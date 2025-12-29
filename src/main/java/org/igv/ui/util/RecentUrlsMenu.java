package org.igv.ui.util;

import org.igv.ui.IGV;
import org.igv.ui.MenuSelectedListener;
import org.igv.ui.RecentUrlsSet;
import org.igv.ui.action.MenuAction;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.List;

public class RecentUrlsMenu extends JMenu {

    public RecentUrlsMenu() {
        this("Recent Files");
    }

    private RecentUrlsMenu(String name) {

        super(name);

        this.addMenuListener((MenuSelectedListener) e -> {
            RecentUrlsSet recentFileList = IGV.getInstance().getRecentUrls();
            if (recentFileList.isEmpty()) {
                RecentUrlsMenu.this.setVisible(false);
            } else {
                RecentUrlsMenu.this.setVisible(true);
                // Remove what's in the menu right now
                RecentUrlsMenu.this.removeAll();
                // Create a menu item for each of the timed autosave files and add it to the menu
                for (ResourceLocator resourceLocator : recentFileList) {
                    JMenuItem menuItem = createMenuItem(resourceLocator);
                    add(menuItem);
                }

                addSeparator();
                JMenuItem clearButton = MenuAndToolbarUtils.createMenuItem(new MenuAction("Clear Recent Files") {
                    @Override
                    public void actionPerformed(ActionEvent event) {
                        IGV.getInstance().getRecentUrls().clear();
                        RecentUrlsMenu.this.setVisible(false);
                    }
                });
                add(clearButton);
            }
        });
    }

    private static JMenuItem createMenuItem(ResourceLocator resourceLocator) {
        MenuAction menuItemAction = new MenuAction(resourceLocator.getPath()) {
            @Override
            public void actionPerformed(ActionEvent event) {
                IGV igv = IGV.getInstance();
                List<ResourceLocator> resource = List.of(resourceLocator);
                igv.loadTracks(resource);
                igv.addToRecentUrls(resource);
            }
        };

        String toolTipText = resourceLocator.getIndexPath() == null
                ? "Load track from " + resourceLocator.getPath()
                : "<html>Load track from <br>path:  " + resourceLocator.getPath()
                    + "<br>index: " + resourceLocator.getIndexPath()
                    + "</html>";

        menuItemAction.setToolTipText(toolTipText);
        return MenuAndToolbarUtils.createMenuItem(menuItemAction);
    }
}

