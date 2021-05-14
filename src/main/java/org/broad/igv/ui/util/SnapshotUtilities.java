/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
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

/**
 * SnapshotUtilities.java
 * <p>
 * Created on November 29, 2007, 2:14 PM
 * <p>
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.log4j.Logger;
import org.broad.igv.ui.panel.Paintable;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

import static org.broad.igv.ui.util.ImageFileTypes.Type.PNG;
import static org.broad.igv.ui.util.ImageFileTypes.Type.SVG;

/**
 * Utility methods for supporting saving of images as jpeg, png, and svg files.
 *
 * @author eflakes
 * @modified jrobinso
 */
public class SnapshotUtilities {

    private static Logger log = Logger.getLogger(SnapshotUtilities.class);

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

        paintImage(target, svgGenerator, width, height, batch);

        Writer out = null;
        try {
            // Finally, stream out SVG to the standard output using
            // UTF-8 encoding.
            boolean useCSS = true; // we want to use CSS style attributes
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

        // Start with a white background
        Color c = g.getColor();
        g.setColor(Color.white);
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
