package org.broad.igv.ui.util;

import org.broad.igv.ui.IGV;
import org.broad.igv.ui.RecentFileSet;
import org.broad.igv.ui.RecentUrlsSet;
import org.broad.igv.ui.action.MenuAction;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import java.awt.event.ActionEvent;
import java.util.List;

public class RecentUrlsMenu extends JMenu {

    public RecentUrlsMenu() {
        this("Recent Files");
    }

    private RecentUrlsMenu(String name) {
        super(name);

        this.addMenuListener(new MenuListener() {
            @Override
            public void menuSelected(MenuEvent e) {
                RecentUrlsSet recentFileList = IGV.getInstance().getRecentUrls();
                if (recentFileList.isEmpty()) {
                    RecentUrlsMenu.this.setVisible(false);
                } else {
                    RecentUrlsMenu.this.setVisible(true);
                    // Remove what's in the menu right now
                    RecentUrlsMenu.this.removeAll();
                    // Create a menu item for each of the timed autosave files and add it to the menu
                    for (ResourceLocator resourceLocator : recentFileList) {
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
                        add(MenuAndToolbarUtils.createMenuItem(menuItemAction));
                    }
                }
            }

            @Override
            public void menuDeselected(MenuEvent e) {
            }

            @Override
            public void menuCanceled(MenuEvent e) {
            }
        });
    }
}

