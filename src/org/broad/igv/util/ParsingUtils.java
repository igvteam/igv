/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.util;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.renderer.*;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * @author jrobinso
 */
public class ParsingUtils {

    private static Logger log = Logger.getLogger(ParsingUtils.class);
    public static final Pattern TAB_PATTERN = Pattern.compile("\t");
    public static final Pattern COMMA_PATTERN = Pattern.compile(",");
    public static final Pattern SEMI_COLON_PATTERN = Pattern.compile(";");
    public static final Pattern EQ_PATTERN = Pattern.compile("=");


    /**
     * Open a BufferedReader on the path, which might be a local file or URL, and might be gzipped or not.
     *
     * @param pathOrUrl
     * @return
     * @throws IOException
     */
    public static BufferedReader openBufferedReader(String pathOrUrl) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(openInputStream(pathOrUrl)));
        return reader;
    }

    public static BufferedReader openBufferedReader(ResourceLocator locator) throws IOException {
        InputStream stream = openInputStreamGZ(locator);
        return new BufferedReader(new InputStreamReader(stream));

    }


    public static AsciiLineReader openAsciiReader(ResourceLocator locator) throws IOException {
        InputStream stream = openInputStreamGZ(locator);
        return new AsciiLineReader(stream);

    }

    public static InputStream openInputStream(String path) throws IOException {
        return openInputStreamGZ(new ResourceLocator(path));
    }

    /**
     * Open an InputStream on the resource.  Wrap it in a GZIPInputStream if neccessary.
     *
     * @param locator
     * @return
     * @throws IOException
     */
    public static InputStream openInputStreamGZ(ResourceLocator locator) throws IOException {

        if (locator.getServerURL() != null) {
            // Use IGV webservice to fetch content
            URL url = new URL(locator.getServerURL() + "?method=getContents&file=" + locator.getPath());
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);
            // Note -- assumption that url stream is compressed!
            try {
                return new GZIPInputStream(is);
            } catch (Exception ex) {
                log.error("Error with gzip stream", ex);
                throw new RuntimeException(
                        "There was a server error loading file: " + locator.getTrackName() +
                                ". Please report to igv-team@broadinstitute.org");

            }

        } else {

            InputStream inputStream = null;
            if (HttpUtils.isRemoteURL(locator.getPath())) {
                URL url = new URL(locator.getPath());
                inputStream = HttpUtils.getInstance().openConnectionStream(url);
            } else {
                String path = locator.getPath();
                if (path.startsWith("file://")) {
                    path = path.substring(7);
                }
                File file = new File(path);
                inputStream = new FileInputStream(file);
            }

            if (locator.getPath().endsWith("gz")) {
                return new GZIPInputStream(inputStream);
            } else {
                return inputStream;
            }
        }
    }

    /**
     * Parse the string and return the result as an integer.  This method supports scientific notation for integers,
     * which Integer.parseInt() does not.
     *
     * @param string
     * @return
     */
    public static int parseInt(String string) {
        return (int) Double.parseDouble(string);
    }


    public static long getContentLength(String path) {
        try {
            long contentLength = -1;
            if (path.startsWith("http:") || path.startsWith("https:")) {
                URL url = new URL(path);
                contentLength = HttpUtils.getInstance().getContentLength(url);

            } else if (path.startsWith("ftp:")) {
                // Use JDK url
                URL url = new URL(path);
                URLConnection connection = url.openConnection();
                connection.setConnectTimeout(Globals.CONNECT_TIMEOUT);
                //For reasons beyond my ken, on Java 7 getContentLength
                //returns -1 without attempting a connection
                //contentLength = connection.getContentLength();
                contentLength = connection.getInputStream().available();
            } else {
                contentLength = (new File(path)).length();
            }
            return contentLength;
        } catch (IOException e) {
            log.error("Error getting content length for: " + path, e);
            return -1;
        }
    }

    public static int estimateLineCount(String path) {

        AsciiLineReader reader = null;
        try {
            final int defaultLength = 100000;
            long fileLength = getContentLength(path);
            if (fileLength <= 0) {
                return defaultLength;
            }

            reader = openAsciiReader(new ResourceLocator(path));
            String nextLine;
            int lines = 0;
            // Skip the first 10 lines (headers, etc)
            int nSkip = 10;
            while (nSkip-- > 0 && reader.readLine() != null) {
            }
            long startPos = reader.getPosition();

            while ((nextLine = reader.readLine()) != null & lines < 100) {
                lines++;
            }

            if (lines == 0) {
                return defaultLength;
            }

            double bytesPerLine = (double) ((reader.getPosition() - startPos) / lines);
            int nLines = (int) (fileLength / bytesPerLine);
            return nLines;

        } catch (Exception e) {
            log.error("Error estimating line count", e);
            return 1000;
        } finally {
            try {
                reader.close();
            } catch (Exception e) {
                // Ignore errors closing reader
            }
        }

    }

    /**
     * Method description
     *
     * @param file
     * @return
     */
    public static List<String> loadRegions(File file) {
        try {
            FileInputStream fileInput = new FileInputStream(file);
            BufferedReader reader = new BufferedReader(new InputStreamReader(fileInput));
            String nextLine;
            List<String> features = new ArrayList<String>();
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                try {
                    if (nextLine.startsWith("chr")) {
                        String[] tokens = nextLine.split("\t");
                        String region = tokens[0] + ":" + tokens[1] + "-" + tokens[2];
                        features.add(region);
                    }
                } catch (NumberFormatException e) {
                    log.error("Error parsing numer in line: " + nextLine);
                }
            }

            reader.close();
            return features;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }


    /**
     * graphType         bar|points           # default is bar
     * yLineMark         real-value           # default is 0.0
     * yLineOnOff        on|off               # default is off
     * windowingFunction maximum|mean|minimum # default is maximum
     * smoothingWindow   off|[2-16]           # default is off
     *
     * @param nextLine
     * @param trackProperties
     * @throws NumberFormatException
     */
    public static boolean parseTrackLine(String nextLine, TrackProperties trackProperties)
            throws NumberFormatException {

        boolean foundProperties = false;
        try {
            // track type=wiggle_0 name="CSF +" description="CSF +" visibility=full autoScale=off viewLimits=-50:50
            List<String> tokens = StringUtils.breakQuotedString(nextLine, ' ');
            for (String pair : tokens) {
                List<String> kv = StringUtils.breakQuotedString(pair, '=');
                if (kv.size() == 2) {
                    foundProperties = true;
                    String key = kv.get(0).toLowerCase().trim();
                    String value = kv.get(1).replaceAll("\"", "");

                    if (key.equals("coords")) {
                        if (value.equals("0")) {
                            trackProperties.setBaseCoord(TrackProperties.BaseCoord.ZERO);
                        } else if (value.equals("1")) {
                            trackProperties.setBaseCoord(TrackProperties.BaseCoord.ONE);
                        }

                    }
                    if (key.equals("name")) {
                        trackProperties.setName(value);
                        //dhmay adding name check for TopHat junctions files. graphType is also checked.
                        if (value.equals("junctions")) {
                            trackProperties.setRendererClass(SpliceJunctionRenderer.class);
                            trackProperties.setHeight(60);
                        }
                    } else if (key.equals("description")) {
                        trackProperties.setDescription(value);
                    } else {
                        final String valueLowerCase = value.toLowerCase();
                        if (key.equals("itemrgb")) {
                            trackProperties.setItemRGB(valueLowerCase.equals("on") || value.equals("1"));
                        } else if (key.equals("usescore")) {
                            trackProperties.setUseScore(value.equals("1"));
                        } else if (key.equals("color")) {
                            Color color = ColorUtilities.stringToColor(value);
                            trackProperties.setColor(color);
                        } else if (key.equals("altcolor")) {
                            Color color = ColorUtilities.stringToColor(value);
                            trackProperties.setAltColor(color);
                        } else if (key.equals("midcolor")) {
                            Color color = ColorUtilities.stringToColor(value);
                            trackProperties.setMidColor(color);
                        } else if (key.equals("autoscale")) {
                            boolean autoscale = value.equals("on");
                            trackProperties.setAutoScale(autoscale);
                        } else if (key.equals("maxheightpixels")) {
                            // There should be 3 values per UCSC spec,  max:default:min.  In the past we have accepted
                            // 2 values,  def:min,  so keep this for backwards compatibility.   IGV currently doesn't
                            // have a "max height"   UPDATE -- apparently 1 value is also allowed!
                            String[] maxDefMin = value.split(":");
                            if (maxDefMin.length >= 2) {
                                int defIDX = (maxDefMin.length == 2 ? 0 : 1);
                                trackProperties.setHeight(Integer.parseInt(maxDefMin[defIDX].trim()));
                                trackProperties.setMinHeight(Integer.parseInt(maxDefMin[defIDX + 1].trim()));
                            } else {
                                // Single value
                                trackProperties.setHeight(Integer.parseInt(value));
                            }

                        } else if (key.equals("url")) {
                            trackProperties.setUrl(value);
                        } else if (key.equals("graphtype")) {

                            if (value.equals("bar")) {
                                trackProperties.setRendererClass(BarChartRenderer.class);
                            } else if (value.equals("points")) {
                                trackProperties.setRendererClass(PointsRenderer.class);
                                trackProperties.setWindowingFunction(WindowFunction.none);
                            } else if (value.equals("line")) {
                                trackProperties.setRendererClass(LineplotRenderer.class);
                            } else if (value.equals("heatmap")) {
                                trackProperties.setRendererClass(HeatmapRenderer.class);
                            } else if (value.equals("junctions")) {
                                //dhmay adding check for graphType=junctions.  name is also checked
                                trackProperties.setRendererClass(SpliceJunctionRenderer.class);
                            } else if (value.equals("genotype")) {
                                trackProperties.setRendererClass(GenotypeRenderer.class);
                            }
                        } else if (key.toLowerCase().equals("viewlimits")) {
                            String[] limits = value.split(":");
                            if (limits.length == 2) {
                                try {
                                    float min = Float.parseFloat(limits[0].trim());
                                    float max = Float.parseFloat(limits[1].trim());
                                    trackProperties.setMinValue(min);
                                    trackProperties.setMaxValue(max);
                                } catch (NumberFormatException e) {
                                    log.error("viewLimits values must be numeric: " + value);
                                }
                            }
                        } else if (key.equals("midrange")) {
                            String[] limits = value.split(":");
                            if (limits.length == 2) {
                                try {
                                    float from = Float.parseFloat(limits[0].trim());
                                    float to = Float.parseFloat(limits[1].trim());
                                    trackProperties.setNeutralFromValue(from);
                                    trackProperties.setNeutralToValue(to);
                                } catch (NumberFormatException e) {
                                    log.error("midrange values must be numeric: " + value);
                                }
                            }
                        } else if (key.equals("ylinemark")) {
                            try {
                                float yLine = Float.parseFloat(value);
                                trackProperties.setyLine(yLine);
                            } catch (NumberFormatException e) {
                                log.error("Number format exception in track line (ylinemark): " + nextLine);
                            }
                        } else if (key.equals("ylineonoff")) {
                            trackProperties.setDrawYLine(value.equals("on"));
                        } else if (key.equals("windowingfunction")) {
                            WindowFunction wf = WindowFunction.getWindowFunction(value);
                            trackProperties.setWindowingFunction(wf);
                        } else if (key.equals("maxfeaturewindow") || key.equals("featurevisibilitywindow") ||
                                key.equals("visibilitywindow")) {
                            try {
                                int windowSize = Integer.parseInt(value);
                                trackProperties.setFeatureVisibilityWindow(windowSize);
                            } catch (NumberFormatException e) {
                                log.error(key + " must be numeric: " + nextLine);

                            }

                        } else if (key.equals("scaletype")) {
                            if (value.equals("log")) {
                                trackProperties.setLogScale(true);
                            }
                        } else if (key.equals("gfftags")) {
                            // Any value other than 0 or off => on
                            boolean gffTags = !(value.equals("0") || (valueLowerCase.equals("off")));
                            trackProperties.setGffTags(gffTags);
                        } else if (key.equals("sortable")) {
                            // Any value other than 0 or off => on
                            boolean sortable = (value.equals("1") || (valueLowerCase.equals("true")));
                            trackProperties.setSortable(sortable);
                        } else if (key.equals("alternateexoncolor")) {
                            trackProperties.setAlternateExonColor(valueLowerCase.equals("on") || value.equals("1"));
                        } else if (key.equals("visibility")) {
                            if (valueLowerCase.equals("1") || valueLowerCase.equals("dense")) {
                                trackProperties.setDisplayMode(Track.DisplayMode.COLLAPSED);
                            } else if (valueLowerCase.equals("2") || valueLowerCase.equals("3") || valueLowerCase.equals("pack")) {
                                trackProperties.setDisplayMode(Track.DisplayMode.EXPANDED);
                            } else if (valueLowerCase.equals("4") || valueLowerCase.equals("squish")) {
                                trackProperties.setDisplayMode(Track.DisplayMode.SQUISHED);
                            }
                        } else if (key.equals("genome") || key.equals("db")) {
                            trackProperties.setGenome(value);
                        } else if (key.equals("bigdataurl") || key.equals("dataurl")) {
                            trackProperties.setDataURL(value);
                        } else if (key.equals("meta")) {
                            trackProperties.setMetaData(value);
                        }
                    }
                }
            }

        } catch (Exception exception) {
            MessageUtils.showMessage("Error parsing track line: " + nextLine + " (" + exception.getMessage() + ")");
        }

        return foundProperties;

    }


    public static boolean pathExists(String covPath) {
        try {
            return (new File(covPath)).exists() ||
                    (HttpUtils.isRemoteURL(covPath) && HttpUtils.getInstance().resourceAvailable(new URL(covPath)));
        } catch (MalformedURLException e) {
            // todo -- log
            return false;
        }
    }
}
