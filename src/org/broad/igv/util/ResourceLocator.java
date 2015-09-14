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

import org.apache.log4j.Logger;
import org.broad.igv.ga4gh.Ga4ghAPIHelper;
import org.broad.igv.gs.GSUtils;
import htsjdk.tribble.Tribble;

import java.awt.*;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Represents a data file or other resource, which might be local file or remote resource.
 *
 * @author jrobinso
 */
public class ResourceLocator {

    private static Logger log = Logger.getLogger(ResourceLocator.class);

    /**
     * Display name
     */
    String name;

    /**
     * The local path or url (http, https, or ftp) for the resource.
     */
    String path;

    /**
     * URL to a database server
     */
    String dbURL;

    /**
     * Optional path to an associated index file
     */
    String indexPath;

    /**
     * URL to a page with general information about the resource
     */
    String trackInforURL;

    /**
     * A URL pattern (UCSC convention) to a specific URL applicable to each feature
     */
    String featureInfoURL;

    /**
     * Descriptive text
     */
    String description;

    /**
     * The type of resource (generally this refers to the file format)
     */
    String type;

    /**
     * Path to an associated density file.  This is used primarily for sequence alignments
     */
    String coverage;

    /**
     * A UCSC style track line.  Overrides value in file, if any.
     */
    String trackLine;  //

    /**
     * Color for features or data.  Somewhat redundant with trackLine.
     */
    Color color;

    String sampleId;

    String username;

    String password;

    private HashMap attributes = new HashMap();

    /**
     * Constructor for local files
     *
     * @param path
     */
    public ResourceLocator(String path) {
        this.setPath(path);
    }

    /**
     * Constructor for database resources
     *
     * @param dbURL
     * @param path
     */
    public ResourceLocator(String dbURL, String path) {
        this.dbURL = dbURL;
        this.setPath(path);
    }

    /**
     * Determines if the resource actually exists.
     *
     * @return true if resource was found.
     */
    public boolean exists() {
        return ParsingUtils.pathExists(path);
    }


    public void setType(String type) {
        this.type = type;
    }

    public String getType() {
        return type;
    }

    /**
     * Return a string suitable for determining file type based on extension
     * May or may not be a full, readable path. txt and gz extensions are stripped
     *
     * @return
     */
    public String getTypeString() {
        if (type != null) {
            return type;
        } else {

            String typeString = path.toLowerCase();
            if (path.startsWith("http://") || path.startsWith("https://")) {
                try {
                    URL url = new URL(path);

                    typeString = url.getPath().toLowerCase();
                    String query = url.getQuery();
                    if (query != null) {
                        Map<String, String> queryMap = HttpUtils.parseQueryString(query);
                        // If type is set explicitly use it
                        if (queryMap.containsKey("dataformat")) {
                            String format = queryMap.get("dataformat");
                            if (format.contains("genomespace")) {
                                typeString = GSUtils.parseDataFormatString(format);
                            } else {
                                typeString = format;
                            }
                        } else if (queryMap.containsKey("file")) {
                            typeString = queryMap.get("file");
                        }
                    }

                } catch (MalformedURLException e) {
                    log.error("Error interpreting url: " + path, e);
                    typeString = path;
                }
            }

            // Strip .txt, .gz, and .xls extensions.  (So  foo.cn.gz => a .cn file)
            if ((typeString.endsWith(".txt") || typeString.endsWith(
                    ".xls") || typeString.endsWith(".gz") || typeString.endsWith(".bgz"))) {
                typeString = typeString.substring(0, typeString.lastIndexOf(".")).trim();
            }

            return typeString;

        }
    }

    /**
     * Returns the portion of the contained path before the query string.
     * If there is no query string, or if the path is not a url,
     * this will be the same as #getPath()
     *
     * @return
     */
    public String getURLPath() {
        return getPath().split("\\?", 2)[0];
    }

    /**
     * Returns the portion of the contained path after the query string.
     * If there is no query string, this will return an empty string
     *
     * @return
     */
    public String getURLQueryString() {
        String[] tmp = getPath().split("\\?", 2);
        if (tmp.length == 1) {
            return "";
        }
        return tmp[1];
    }

    public String toString() {
        return path + (dbURL == null ? "" : " " + dbURL);
    }

    public String getPath() {
        return path;
    }

    public String getFileName() {
        return (new File(path)).getName();
    }


    public String getDBUrl() {
        return dbURL;
    }

    public boolean isLocal() {
        return dbURL == null && !FileUtils.isRemote(path) && !Ga4ghAPIHelper.RESOURCE_TYPE.equals(type);
    }

    public void setTrackInforURL(String trackInforURL) {
        this.trackInforURL = trackInforURL;
    }

    public String getTrackInfoURL() {
        return trackInforURL;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getTrackName() {
        if(name == null) {
            if(path.startsWith("http://") || path.startsWith("https://")) {
                try {
                    return new File((new URL(path)).getPath()).getName();
                } catch (MalformedURLException e) {
                    return path;
                }
            }
            else {
                return new File(path).getName();
            }
        }
        return name;

    }


    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public String getCoverage() {
        return coverage;
    }

    public void setCoverage(String coverage) {
        this.coverage = coverage;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }


    public String getFeatureInfoURL() {
        return featureInfoURL;
    }

    public void setFeatureInfoURL(String featureInfoURL) {
        this.featureInfoURL = featureInfoURL;
    }

    public void setPath(String path) {
        if (path != null && path.startsWith("file://")) {
            this.path = path.substring(7);
        } else {
            this.path = path;
        }
    }

    public String getTrackLine() {
        return trackLine;
    }

    public void setTrackLine(String trackLine) {
        this.trackLine = trackLine;
    }

    public String getSampleId() {
        return sampleId;
    }

    public void setSampleId(String sampleId) {
        this.sampleId = sampleId;
    }

    public String getIndexPath() {
        return indexPath;
    }

    public void setIndexPath(String indexPath) {
        this.indexPath = indexPath;
    }

    public String getUsername() {
        return username;
    }

    public void setUsername(String username) {
        this.username = username;
    }

    public String getPassword() {
        return password;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ResourceLocator that = (ResourceLocator) o;

        if (dbURL != null ? !dbURL.equals(that.dbURL) : that.dbURL != null) return false;
        if (path != null ? !path.equals(that.path) : that.path != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = path != null ? path.hashCode() : 0;
        result = 31 * result + (dbURL != null ? dbURL.hashCode() : 0);
        return result;
    }

    public String getBamIndexPath() {

        if (indexPath != null) return indexPath;

        if (path.toLowerCase().startsWith("http://") || path.toLowerCase().startsWith("https://")) {
            // See if bam file is specified by parameter
            try {
                URL url = new URL(path);
                String queryString = url.getQuery();
                if (queryString != null) {
                    Map<String, String> parameters = HttpUtils.parseQueryString(queryString);
                    if (parameters.containsKey("index")) {
                        return parameters.get("index");
                    } else if (parameters.containsKey("file")) {
                        String bamFile = parameters.get("file");
                        String bamIndexFile = bamFile + ".bai";
                        String newQueryString = queryString.replace(bamFile, bamIndexFile);
                        return path.replace(queryString, newQueryString);
                    }
                    else {
                        String ip = path.replace(url.getPath(), url.getPath() + ".bai");
                        return ip;
                    }
                }
            } catch (MalformedURLException e) {
                log.error(e.getMessage(), e);
            }

        }

        return path + ".bai";
    }

    /**
     * Add the {@code indexExtension} to the path in locator, preserving
     * query string elements if present
     *
     * @param locator
     * @param indexExtension
     * @return
     */
    public static String appendToPath(ResourceLocator locator, String indexExtension) {
        String indexFile = locator.getURLPath() + indexExtension;
        String qs = locator.getURLQueryString();
        if (qs != null && qs.length() > 0) {
            indexFile += "?" + qs;
        }
        return indexFile;
    }

    /**
     * @param locator
     * @return locator.getIndexPath() if not null, otherwise
     * {@link #appendToPath(ResourceLocator, String)}
     * where the second argument is .idx or tbi, depending on the resource
     */
    public static String indexFile(ResourceLocator locator) {
        if (locator.getIndexPath() != null) {
            return locator.getIndexPath();
        }
        String indexExtension = (locator.getPath().toLowerCase().endsWith(".gz") || locator.getPath().toLowerCase().endsWith(".bgz")) ? ".tbi" : Tribble.STANDARD_INDEX_EXTENSION;
        return appendToPath(locator, indexExtension);
    }

    public void setAttribute(String key, Object value) {
        this.attributes.put(key, value);
    }

    public Object getAttribute(String key) {
        return attributes.get(key);
    }


    /**
     * FOR LOAD FROM SERVER
     */
    public static enum AttributeType {

        DB_URL("serverURL"),
        PATH("path"),
        DESCRIPTION("description"),
        HYPERLINK("hyperlink"),
        INFOLINK("infolink"),
        ID("id"),
        SAMPLE_ID("sampleId"),
        NAME("name"),
        URL("url"),
        RESOURCE_TYPE("resourceType"),
        TRACK_LINE("trackLine"),
        COVERAGE("coverage"),
        COLOR("color");

        private String name;

        AttributeType(String name) {
            this.name = name;
        }

        public String getText() {
            return name;
        }

        @Override
        public String toString() {
            return getText();
        }

    }
}
