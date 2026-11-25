package org.igv.ui;

import java.awt.Desktop;
import java.awt.Image;
import java.awt.Taskbar;

/**
  * Java version-specific integration with the platform Desktop and particularly 
  * for OS X (macOS) specific items. * @author eby
 */
public class DesktopIntegration {
    public static void setDockIcon(Image image) {
        Taskbar.getTaskbar().setIconImage(image);
    }

    public static void setAboutHandler(IGVMenuBar igvMenuBar) {
        Desktop.getDesktop().setAboutHandler(e -> igvMenuBar.showAboutDialog());
    }

    public static void setQuitHandler() {
        Desktop.getDesktop().setQuitHandler((e, response) -> {
            response.performQuit();
        });
    }
}
