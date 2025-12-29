package org.broad.igv.ui;

import java.awt.Desktop;
import java.awt.GraphicsEnvironment;
import java.awt.Image;
import java.awt.Taskbar;

import javax.swing.JOptionPane;

/**
  * Java version-specific integration with the platform Desktop and particularly 
  * for OS X (macOS) specific items. * @author eby
 */
public class DesktopIntegration {
    public static final void verifyJavaPlatform() {
        String javaVersion = System.getProperty("java.version");
        if (javaVersion == null || javaVersion.startsWith("1.8")) {
            try {
                System.out.println("Detected an unsupported Java version.  Java 8 is not supported by this release.");

                if (!GraphicsEnvironment.isHeadless()) {
                    JOptionPane.showMessageDialog(null, "Detected an unsupported Java version.  Java 8 is not supported by this release.");
                }
            } finally {
                System.exit(1);
            }
        }
    }
    
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
