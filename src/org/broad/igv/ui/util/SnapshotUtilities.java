/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/**
 * SnapshotUtilities.java
 *
 * Created on November 29, 2007, 2:14 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.Paintable;
import org.broad.igv.ui.svg.SVGGraphics;
import org.broad.igv.util.Utilities;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.*;
import java.util.LinkedHashMap;

/**
 * Utility methods for supporting saving of images as jpeg, png, and svg files.
 *
 * @author eflakes
 * @modified jrobinso
 */
public class SnapshotUtilities {

    /**
     * Class logger
     */
    private static Logger logger = Logger.getLogger(SnapshotUtilities.class);

    private static LinkedHashMap<SnapshotFileType, SnapshotFileFilter> SNAPSHOT_TYPE_TO_FILTER = new LinkedHashMap();

    /**
     * The maximum height in pixels for snapshots of a panel.
     */
    public static int MAX_PANEL_HEIGHT = 1000;


    static {

        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.JPEG, new SnapshotFileFilter(SnapshotFileType.JPEG));
        //SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.PDF,
        //    new SnapshotFileFilter(SnapshotFileType.PDF));
        //SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.EPS,
        //        new SnapshotFileFilter(SnapshotFileType.EPS));
        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.SVG, new SnapshotFileFilter(SnapshotFileType.SVG));
        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.PNG, new SnapshotFileFilter(SnapshotFileType.PNG));
    }

    public static FileFilter[] getAllSnapshotFileFilters() {
        return SNAPSHOT_TYPE_TO_FILTER.values().toArray(new FileFilter[SNAPSHOT_TYPE_TO_FILTER.size()]);
    }

    public static SnapshotFileFilter getSnapshotFileFilterForType(SnapshotFileType type) {
        return SNAPSHOT_TYPE_TO_FILTER.get(type);
    }

    /**
     * Creates a device compatible buffered svg.
     *
     * @param width  the svg width in pixels
     * @param height the svg height in pixels
     */
    public static BufferedImage getDeviceCompatibleImage(int width, int height) {

        GraphicsEnvironment graphicsEnvironment = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice screenDevice = graphicsEnvironment.getDefaultScreenDevice();
        GraphicsConfiguration graphicConfiguration = screenDevice.getDefaultConfiguration();
        BufferedImage image = graphicConfiguration.createCompatibleImage(width, height);

        return image;
    }

    /**
     * Snapshot types
     */
    public static enum SnapshotFileType {

        NULL("", ""),
        EPS(".eps", "Encapsulated Postscript Files (*.eps)"),
        PDF(".pdf", "Portable Document FormatFles (*.pdf)"),
        SVG(".svg", "Scalable Vector Graphics Files (*.svg)"),
        PNG(".png", "Portable Network Graphics Files (*.png)"),
        JPEG(".jpeg", "Joint Photographic Experts Group Files (*.jpeg)");
        private String fileExtension;
        private String fileDescription;

        SnapshotFileType(String extension, String description) {
            fileExtension = extension;
            fileDescription = description;
        }

        public String getExtension() {
            return fileExtension;
        }

        public String getDescription() {
            return fileDescription;
        }
    }

    public static String getFileExtension(String filePath) {

        String extension = null;

        int indexOfExtension = filePath.lastIndexOf(".");
        if (indexOfExtension >= 0) {
            extension = filePath.substring(indexOfExtension, filePath.length());
        }
        return extension;
    }

    public static SnapshotFileType getSnapshotFileType(String fileExtension) {

        String extension = fileExtension.toLowerCase();
        SnapshotFileType type = null;

        if (SnapshotFileType.EPS.getExtension().equals(extension)) {
            type = SnapshotFileType.EPS;
        } else if (SnapshotFileType.PDF.getExtension().equals(extension)) {
            type = SnapshotFileType.PDF;
        } else if (SnapshotFileType.SVG.getExtension().equals(extension)) {
            type = SnapshotFileType.SVG;
        } else if (SnapshotFileType.PNG.getExtension().equals(extension)) {
            type = SnapshotFileType.PNG;
        } else if (SnapshotFileType.JPEG.getExtension().equals(extension)) {
            type = SnapshotFileType.JPEG;
        } else {
            type = SnapshotFileType.NULL;
        }
        return type;
    }


    /**
     * Snapshot file filter
     */
    public static class SnapshotFileFilter extends FileFilter {

        private SnapshotFileType type = SnapshotFileType.EPS;

        public SnapshotFileFilter(SnapshotFileType type) {
            this.type = type;
        }

        public boolean accept(File file) {

            if (file.isDirectory()) {
                return true;
            }

            return file.getName().toLowerCase().endsWith(type.getExtension());
        }

        public String getDescription() {
            return type.getDescription();
        }

        public String getExtension() {
            return type.getExtension();
        }

        public boolean accept(File file, String name) {
            return name.toLowerCase().endsWith(type.getExtension());
        }
    }

    public static void doComponentSnapshot(Component component, File file, SnapshotFileType type) {

        int width = component.getWidth();
        int height = component.getHeight();

        // Call appropriate converter
        switch (type) {
            case JPEG:
                exportScreenShotJPEG(component, file, width, height);
                break;
            //case EPS:
            //    exportScreenShotEPS(component, file, width, height);
            //    break;
            case PNG:
                exportScreenShotPNG(component, file, width, height);
                break;
            case SVG:
                logger.debug("Exporting svg screenshot");
                exportScreenshotSVG(component, file);
                break;
        }
    }

    private static void exportScreenshotSVG2(Component target, File selecteddFile) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(selecteddFile));

            pw.println("<?xml version=\"1.0\" standalone=\"no\"?>\n" +
                    "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n" +
                    "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" +
                    "\n" +
                    "<svg width=\"100%\" height=\"100%\" version=\"1.1\"\n" +
                    "xmlns=\"http://www.w3.org/2000/svg\">");


            SVGGraphics g2d = new SVGGraphics(pw);

            target.paint(g2d);
            pw.print("</svg>");
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }


    private static void exportScreenshotSVG(Component target, File selectedFile) {
        // Get a DOMImplementation.
        try {
            logger.debug("Getting dom");
            DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();


            // Create an instance of org.w3c.dom.Document.
            String svgNS = "http://www.w3.org/2000/svg";
            Document document = domImpl.createDocument(svgNS, "svg", null);


            // Create an instance of the SVG Generator.
            SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
            //logger.info("Painting");
            target.paintAll(svgGenerator);

            // Finally, stream out SVG to the standard output using
            // UTF-8 encoding.
            boolean useCSS = true; // we want to use CSS style attributes
            Writer out = new BufferedWriter(new FileWriter(selectedFile));
            //logger.info("Writing output");
            svgGenerator.stream(out, useCSS);
            //logger.info("Done");
        } catch (Exception e) {
            logger.error("Error creating SVG file", e);
            MessageUtils.showMessage("Error encountered creating SVG file: " + e.toString());
        }
    }

    private static void exportScreenShotJPEG(Component target, File selectedFile, int width, int height) {

        BufferedImage image = getDeviceCompatibleImage(width, height); //  new BufferedImage(width, height, BufferedImage.TYPE_BYTE_INDEXED);
        Graphics g = image.createGraphics();
        target.paintAll(g);

        if (selectedFile != null) {

            if (!selectedFile.getName().toLowerCase().endsWith(".jpeg")) {
                String correctedFilename = selectedFile.getAbsolutePath() + ".jpeg";
                selectedFile = new File(correctedFilename);
            }
            writeImage(image, selectedFile, "jpeg");
        }
    }

    private static void exportScreenShotPNG(Component target, File selectedFile, int width, int height) {

        BufferedImage image = getDeviceCompatibleImage(width, height);
        Graphics g = image.createGraphics();
        target.paintAll(g);

        if (selectedFile != null) {

            if (!selectedFile.getName().toLowerCase().endsWith(".png")) {
                String correctedFilename = selectedFile.getAbsolutePath() + ".png";
                selectedFile = new File(correctedFilename);
            }
            writeImage(image, selectedFile, "png");
        }
    }


    public static String doSnapshotOffscreen(Component target, File selectedFile) {

        if (!(target instanceof Paintable)) {
            // TODO -- message that target does not support this
            return "Error: target is not paintable";
        }

        try {
            IGV.getInstance().setExportingSnapshot(true);

            Rectangle rect = target.getBounds();

            int height = ((MainPanel) target).getOffscreenImageHeight();
            rect.height = height;

            // translate to (0, 0) if neccessary
            int dx = rect.x;
            int dy = rect.y;
            rect.x = 0;
            rect.y = 0;
            rect.width -= dx;
            rect.height -= dy;

            Paintable paintable = (Paintable) target;

            final String filenameLowercase = selectedFile.getName().toLowerCase();
            if (filenameLowercase.endsWith(".svg")) {
                logger.debug("Getting dom");
                DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();


                // Create an instance of org.w3c.dom.Document.
                String svgNS = "http://www.w3.org/2000/svg";
                Document document = domImpl.createDocument(svgNS, "svg", null);

                // Create an instance of the SVG Generator.
                SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
                //logger.info("Painting");
                paintable.paintOffscreen(svgGenerator, rect);

                // Finally, stream out SVG to the standard output using
                // UTF-8 encoding.
                boolean useCSS = true; // we want to use CSS style attributes
                Writer out = new BufferedWriter(new FileWriter(selectedFile));
                //logger.info("Writing output");
                svgGenerator.stream(out, useCSS);


            } else {

                BufferedImage image = getDeviceCompatibleImage(rect.width, rect.height);
                Graphics2D g = image.createGraphics();
                paintable.paintOffscreen(g, rect);
                if (selectedFile != null) {
                    if (!filenameLowercase.endsWith(".png")) {
                        String correctedFilename = selectedFile.getAbsolutePath() + ".png";
                        selectedFile = new File(correctedFilename);
                    }
                    writeImage(image, selectedFile, "png");
                }
            }

        } catch (Exception e) {
            logger.error("Error creating snapshot", e);
            return "Error: " + e.toString();
        } finally {
            IGV.getInstance().setExportingSnapshot(false);
        }
        return "OK";

    }

    private static void writeImage(BufferedImage image, File f, String type) {
        try {
            ImageIO.write(image, type, f);
        } catch (IOException e) {
            logger.error(("Error creating: " + f.getAbsolutePath()), e);
        }
    }


    public static class SnapshotFileChooser extends JFileChooser {

        boolean accepted = false;
        File previousFile;

        public SnapshotFileChooser(File directory, File selectedFile) {
            super(directory);
            setPreviousFile(selectedFile);
            init();
        }


        public void approveSelection() {
            accepted = true;
            super.approveSelection();
        }

        public void setPreviousFile(File file) {
            this.previousFile = file;
            setSelectedFile(previousFile);
        }

        public File getPreviousFile() {
            return previousFile;
        }

        @Override
        public void cancelSelection() {
            setSelectedFile(null);
            super.cancelSelection();
        }


        @Override
        protected JDialog createDialog(Component parent) throws HeadlessException {
            JDialog dialog = super.createDialog(parent);
            dialog.setLocation(300, 200);
            dialog.setResizable(false);
            dialog.addWindowListener(new WindowAdapter() {

                @Override
                public void windowClosing(WindowEvent e) {
                    if (!accepted) {
                        setSelectedFile(null);
                    }
                }
            });
            return dialog;
        }

        private void init() {

            FileFilter[] fileFilters = SnapshotUtilities.getAllSnapshotFileFilters();
            // Setup FileFilters
            if (fileFilters != null) {
                for (FileFilter fileFilter : fileFilters) {
                    addChoosableFileFilter(fileFilter);
                }
            }

            addPropertyChangeListener(new PropertyChangeListener() {
                public void propertyChange(PropertyChangeEvent e) {


                    File oldFile = null;
                    String property = e.getPropertyName();
                    if (JFileChooser.SELECTED_FILE_CHANGED_PROPERTY.equals(property)) {
                        oldFile = (File) e.getOldValue();
                    } else if (JFileChooser.FILE_FILTER_CHANGED_PROPERTY.equals(property)) {

                        if (e.getOldValue() instanceof SnapshotFileFilter &&
                                e.getNewValue() instanceof SnapshotFileFilter) {

                            SnapshotFileFilter newFilter = (SnapshotFileFilter) e.getNewValue();

                            File currentDirectory = getCurrentDirectory();
                            File previousFile = getPreviousFile();
                            if (previousFile != null) {

                                File file = null;
                                if (currentDirectory != null) {
                                    file = new File(currentDirectory, previousFile.getName());
                                } else {
                                    file = previousFile;
                                }

                                final File selectedFile = Utilities.changeFileExtension(
                                        file, newFilter.getExtension());

                                UIUtilities.invokeOnEventThread(new Runnable() {

                                    public void run() {
                                        setPreviousFile(selectedFile);
                                        validate();
                                    }
                                });
                            }

                        }
                    }
                }
            });
        }
    }
}
