/**
 * SnapshotUtilities.java
 * <p>
 * Created on November 29, 2007, 2:14 PM
 * <p>
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.igv.ui.util;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.igv.logging.*;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.Paintable;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

import static org.igv.ui.util.ImageFileTypes.Type.PNG;
import static org.igv.ui.util.ImageFileTypes.Type.SVG;

/**
 * Utility methods for supporting saving of images as jpeg, png, and svg files.
 *
 * @author eflakes
 * @modified jrobinso
 */
public class SnapshotUtilities {

    private static Logger log = LogManager.getLogger(SnapshotUtilities.class);

    /**
     * The maximum height in pixels for snapshots of a panel.
     */
    public static int DEFAULT_MAX_PANEL_HEIGHT = 1000;

    /**
     * We need to use a static for max panel height,  or alternatively much refactoring
     */
    private static int maxPanelHeight = DEFAULT_MAX_PANEL_HEIGHT;

    public static int getMaxPanelHeight() {
        return maxPanelHeight;
    }

    public static void setMaxPanelHeight(int h) {
        maxPanelHeight = h;
    }

    public static void resetMaxPanelHeight() {
        maxPanelHeight = DEFAULT_MAX_PANEL_HEIGHT;
    }

    public static boolean snapshotInProgress = false;

    // Treat this class as a singleton, no instances allowed
    private SnapshotUtilities() {
    }


    public static String doComponentSnapshot(Component component, File file, ImageFileTypes.Type type, boolean batch) throws IOException {

        try {
            snapshotInProgress = true;
            if (!(component instanceof Paintable)) {
                throw new RuntimeException("Error: " + component + " is not an instance of Paintable");
            }

            Paintable paintable = (Paintable) component;
            int width = component.getWidth();
            int height = paintable.getSnapshotHeight(batch);

            // Call appropriate converter
            if (type == SVG) {
                exportScreenshotSVG((Paintable) component, file, width, height, batch);
                return "OK";
            } else if (type == PNG) {
                String format = "png";
                String[] exts = new String[]{"." + format};
                exportScreenShotBufferedImage((Paintable) component, file, width, height, exts, format, batch);
                return "OK";
            } else {
                final String message = "No image write for file type: " + file + " Try '.png' or '.svg'";
                MessageUtils.showMessage(message);
                return "ERROR: " + message;
            }
        } finally {
            snapshotInProgress = false;
        }
    }

    private static void exportScreenshotSVG(Paintable target, File selectedFile, int width, int height, boolean batch) throws IOException {

        String format = "svg";
        selectedFile = fixFileExt(selectedFile, new String[]{format}, format);

        // Create an instance of org.w3c.dom.Document.                                                                                      
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, format, null);

        // Write image data into document                                                                                                   
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        final Color background = UIConstants.getTrackPanelBackground();
        svgGenerator.setBackground(background);
        svgGenerator.setColor(background);
        svgGenerator.fillRect(0, 0, width, height);

        paintImage(target, svgGenerator, width, height, batch);

        Writer out = null;
        try {
            // Finally, stream out SVG to the standard output using                                                                         
            // UTF-8 encoding.                                                                                                              
            boolean useCSS = false; // we want to use CSS style attributes
            out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(selectedFile), "UTF-8"));
            svgGenerator.stream(out, useCSS);
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                log.error("Error closing svg file", e);
            }
        }
    }

    /**
     * Export the specified {@code target} component as a {@code BufferedImage} to the given file.
     *
     * @param target
     * @param selectedFile
     * @param width
     * @param height
     * @param allowedExts
     * @param format       Format, also appended as an extension if the file doesn't end with anything in {@code allowedExts}
     * @throws IOException
     */
    private static void exportScreenShotBufferedImage(Paintable target, File selectedFile, int width, int height,
                                                      String[] allowedExts, String format, boolean batch) throws IOException {

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = image.createGraphics();

        // Start with the panel background rather than a hardcoded white -- paintOffscreen implementations do not
        // fill a background of their own, so in dark mode the export was a white page under white text.  Setting
        // the graphics background as well as filling lets renderers that erase behind text match it.
        final Color background = UIConstants.getTrackPanelBackground();
        Color c = g.getColor();
        g.setBackground(background);
        g.setColor(background);
        g.fillRect(0, 0, width, height);
        g.setColor(c);

        paintImage(target, g, width, height, batch);

        selectedFile = fixFileExt(selectedFile, allowedExts, format);
        if (selectedFile != null) {
            log.debug("Writing image to " + selectedFile.getAbsolutePath());
            boolean success = ImageIO.write(image, format, selectedFile);
            if (!success) {
                MessageUtils.showMessage("Error writing image file of type: " + format + ". Try .png or .svg");
            }
        }
    }


    private static void paintImage(Paintable target, Graphics2D g, int width, int height, boolean batch) {
        Rectangle rect = new Rectangle(0, 0, width, height);
        target.paintOffscreen(g, rect, batch);
    }


    /**
     * Add a file extension to the file if it doesn't already
     * have an acceptable one
     *
     * @param selectedFile
     * @param allowedExts  Strings which qualify as extensions
     * @param defExtension Default extension. A period be inserted in between the file path iff {@code defExtension}
     *                     does not already have it
     * @return Either the input File, if it had an extension contained in {@code allowedExts},
     * or a new with with {@code defExtension} appended
     */
    private static File fixFileExt(File selectedFile, String[] allowedExts, String defExtension) {
        boolean hasExt = false;
        if (selectedFile != null) {
            for (String ext : allowedExts) {
                if (selectedFile.getName().toLowerCase().endsWith(ext)) {
                    hasExt = true;
                    break;
                }
            }
            if (!hasExt) {
                String addExt = defExtension.startsWith(".") ? defExtension : "." + defExtension;
                String correctedFilename = selectedFile.getAbsolutePath() + addExt;
                selectedFile = new File(correctedFilename);
            }
        }
        return selectedFile;
    }


    /**
     * Paint a live component hierarchy into a PNG using Swing's ordinary paint path.
     * <p>
     * This is a debugging aid, and differs from the image export above in two ways that matter when reviewing
     * the UI itself rather than the tracks:
     * <ul>
     *   <li>It renders whatever is actually on screen -- scrollbars, checkboxes, the track selection strip,
     *       borders, dialogs -- rather than only the components implementing {@link Paintable}, whose
     *       {@code paintOffscreen} methods deliberately skip the surrounding chrome.</li>
     *   <li>It touches no screen capture API, so it needs no macOS "Screen Recording" permission and works
     *       over a remote or automated session.</li>
     * </ul>
     * What it cannot capture is anything the OS draws rather than Swing: the macOS screen menu bar, native
     * window chrome, and heavyweight popups.
     *
     * @param component the component to paint, typically {@code IGV.getInstance().getContentPane()}
     * @param file      destination PNG
     * @param scale     integer scale factor; use 2 for legible text when inspecting colors
     */
    public static void writeComponentImage(Component component, File file, int scale) throws IOException {

        if (component.getWidth() <= 0 || component.getHeight() <= 0) {
            throw new IOException("Component has not been laid out: " + component.getClass().getSimpleName());
        }
        if (scale < 1) {
            scale = 1;
        }

        final int s = scale;
        final BufferedImage image = new BufferedImage(component.getWidth() * s, component.getHeight() * s,
                BufferedImage.TYPE_INT_RGB);

        // Swing components may only be painted on the event dispatch thread.  printAll() rather than paint() so
        // double buffering is bypassed -- painting a buffered component directly can yield a blank image.
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            Graphics2D g = image.createGraphics();
            try {
                g.scale(s, s);
                g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                component.printAll(g);
            } finally {
                g.dispose();
            }
        });

        File parent = file.getAbsoluteFile().getParentFile();
        if (parent != null && !parent.exists()) {
            parent.mkdirs();
        }
        ImageIO.write(image, "png", file);
        log.info("Wrote component image to " + file.getAbsolutePath());
    }

    /**
     * Creates a device compatible BufferedImage
     *
     * @param width  the width in pixels
     * @param height the height in pixels
     */
    public static BufferedImage getDeviceCompatibleImage(int width, int height) {

        GraphicsEnvironment graphicsEnvironment = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice screenDevice = graphicsEnvironment.getDefaultScreenDevice();
        GraphicsConfiguration graphicConfiguration = screenDevice.getDefaultConfiguration();
        BufferedImage image = graphicConfiguration.createCompatibleImage(width, height);

        return image;
    }


}
