/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

package org.broad.igv.util;

import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPReply;
import org.broad.igv.util.ftp.FTPUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.renderer.*;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.MessageUtils;
import htsjdk.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
    public static final Pattern PERIOD_PATTERN = Pattern.compile("\\.");
    public static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s+");


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
     * Open an InputStream on the resource.  Wrap it in a GZIPInputStream if necessary.
     *
     * @param locator
     * @return
     * @throws IOException
     */
    public static InputStream openInputStreamGZ(ResourceLocator locator) throws IOException {

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

    /**
     * Load a text file of format
     * {@code key}={@code value}
     * <p/>
     * Lines beginning with a "#" are skipped
     * <p/>
     * The {@code value} field must not contain any "=", no allowance is made for quoting or escaping
     * or anything like that.
     *
     * @param inputStream
     * @return
     */
    public static Map<String, String> loadMap(InputStream inputStream) {
        BufferedReader reader = null;
        Map<String, String> map = new HashMap<String, String>();
        try {
            reader = new BufferedReader(new InputStreamReader(inputStream));
            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {
                if (nextLine.startsWith("#")) continue;

                String[] tokens = nextLine.split("=");
                if (tokens.length == 2) {
                    map.put(tokens[0], tokens[1]);
                } else {
                    throw new IllegalArgumentException("Incorrect number of tokens at line: " + nextLine);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ex) {
                log.error(ex.getMessage(), ex);
            }
        }

        return map;
    }

    private static final String codeFilePath = "resources/iupac_regex_table.txt";

    public static Map<String, String> loadIUPACMap() {
        return loadMap(ParsingUtils.class.getResourceAsStream(codeFilePath));
    }


    private static final DateFormat ftpDateFormat = new SimpleDateFormat("yyyyMMddHHmmss");

    /**
     * Returns the number of milliseconds since January 1, 1970, 00:00:00 GMT
     * since the specified resource was modified, or 0 if not known/error
     *
     * @param path
     * @return
     */
    public static long getLastModified(String path) {
        if (HttpUtils.isRemoteURL(path)) {
            String resp = null;
            try {
                URL url = new URL(path);
                if (path.startsWith("ftp:")) {
                    String host = url.getHost();
                    FTPClient ftp = FTPUtils.connect(host, url.getUserInfo(), new UserPasswordInputImpl());
                    ftp.pasv();
                    FTPReply reply = ftp.executeCommand("MDTM " + url.getPath());
                    resp = reply.getReplyString();
                    return ftpDateFormat.parse(resp).getTime();
                } else {
                    return HttpUtils.getInstance().getLastModified(url);
                }

            } catch (MalformedURLException e) {
                log.error("Malformed url " + path, e);
            } catch (IOException e) {
                log.error("Error getting modified date for " + path, e);
            } catch (ParseException e) {
                log.error("Error parsing Last-Modified " + resp, e);
            } catch (NumberFormatException e) {
                log.error("Error parsing Last-Modified " + resp, e);
            }
            return 0;

        } else {
            File f = new File(path);
            return f.exists() ? f.lastModified() : 0;
        }
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
        if (covPath == null) return false;
        try {
            return (new File(covPath)).exists() ||
                    (HttpUtils.isRemoteURL(covPath) && HttpUtils.getInstance().resourceAvailable(new URL(covPath)));
        } catch (MalformedURLException e) {
            log.error(e.getMessage(), e);
            return false;
        }
    }

    /**
     * Return the "IGV extension" (basically the extension after stripping trailing qualifiers) for the input path.
     * his is the string IGV uses to identify the format and data type of the file.
     *
     * @param path
     * @return
     */
    public static String getIGVExtension(String path) {

        // String off gzip first
        if (path.endsWith(".gz")) path = path.substring(0, path.length() - 3);

        // Now common qualifiers
        if (path.endsWith(".txt") || path.endsWith(".xls")) path = path.substring(0, path.length() - 4);

        int idx = path.lastIndexOf('.');
        return idx < 0 ? path : path.substring(idx + 1, path.length());
    }
}
