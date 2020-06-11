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

import com.google.gson.JsonObject;
import htsjdk.tribble.Tribble;
import org.apache.log4j.Logger;
import org.broad.igv.google.Ga4ghAPIHelper;
import org.broad.igv.google.GoogleUtils;

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
     * /**
     * Path to an associated density file.  This is used primarily for sequence alignments
     */
    String coverage;

    /**
     * Optional path to an associated variant->bam mapping file (vcf only)
     */

    String mappingPath;

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
     * A UCSC style track line.  Overrides value in file, if any.
     */
    String trackLine;  //

    /**
     * Color for features or data.  Somewhat redundant with trackLine.
     */
    Color color;

    String sampleId;
    
    private HashMap attributes = new HashMap();
    private boolean indexed;

    /**
     * Constructor for local files
     *
     * @param path
     */
    public ResourceLocator(String path) {
        this.setPath(path);

        if (path != null && path.startsWith("https://") && GoogleUtils.isGoogleDrive(path)) {
            this.resolveGoogleDrive(path);
        }

    }

    private void resolveGoogleDrive(String path) {

        JsonObject fileInfo = GoogleUtils.getDriveFileInfo(path);
        this.name = fileInfo.get("name").getAsString();
        this.type = getTypeString(this.name);
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
        return ParsingUtils.fileExists(path);
    }


    public void setType(String type) {
        this.type = type;
    }

    public String getType() {
        return type;
    }

    public String getTypeString() {
        return this.getTypeString(this.path);
    }

    /**
     * Return a string suitable for determining file type based on extension
     * May or may not be a full, readable path. txt and gz extensions are stripped
     *
     * @return
     */
    public String getTypeString(String pathOrName) {
        if (type != null) {
            return type;
        } else {

            String typeString = pathOrName.toLowerCase();
            if (typeString.startsWith("http://") || typeString.startsWith("https://") ||
                    typeString.startsWith("gs://") || typeString.startsWith("s3://")) {

                try {
                    URL url = HttpUtils.createURL(pathOrName);

                    typeString = url.getPath().toLowerCase();
                    String query = url.getQuery();
                    if (query != null) {
                        Map<String, String> queryMap = URLUtils.parseQueryString(query);
                        // If type is set explicitly use it
                        if (queryMap.containsKey("dataformat")) {
                            String format = queryMap.get("dataformat");
                            typeString = format;
                        } else if (queryMap.containsKey("file")) {
                            typeString = queryMap.get("file");
                        }
                    }

                } catch (MalformedURLException e) {
                    log.error("Error interpreting url: " + pathOrName, e);
                    typeString = pathOrName;
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
        if (name == null) {
            if (path.startsWith("http://") || path.startsWith("https://") || path.startsWith("gs://")) {
                int idx = path.lastIndexOf('/');
                int idx2 = path.indexOf('?');
                return idx2 > idx ? path.substring(idx + 1, idx2) : path.substring(idx + 1);

            } else {
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
            this.path = path.substring("file://".length());
        } else if (path != null && path.startsWith("gs://")) {
            this.path = GoogleUtils.translateGoogleCloudURL(path);
        } else if (path != null && path.startsWith("s3://")) {
            this.path = path;

            // Set UI human-readable short name for the file
            String objFname = "";
            if (path.contains("/")) {
                objFname = path.substring(path.lastIndexOf('/')).replace("/", "");
            } else {
                objFname = path;
            }

            log.debug("S3 object filename visible in IGV UI is: " + objFname);
            this.setName(objFname);

            String s3UrlIndexPath = detectIndexPath(path);

            this.setIndexPath(s3UrlIndexPath);

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

    // XXX: Why does IGV not do that across all providers already?

    /**
     * Takes in a non-pre-signed URL and returns its (guessed) indexfile.
     *
     * @param inputPath: Path containing vcf/bam file
     * @return indexPath: Guessed path containing the corresponding index (in the CWD-equivalent dir level)
     */
    public String detectIndexPath(String inputPath) {
        log.debug("detectIndexPath() input S3 path is: " + inputPath);
        String indexPath = "";
        if (inputPath.contains(".bam")) {
            indexPath = inputPath + ".bai";
        } else if (inputPath.endsWith(".gz")) {
            indexPath = inputPath + ".tbi";
        } else {
            log.debug("S3 index object filetype could not be determined from S3 url");
        }
        return indexPath;
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

    public String getMappingPath() {
        return mappingPath;
    }

    public void setMappingPath(String mappingPath) {
        this.mappingPath = mappingPath;
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
        } else {
            if (isCloudOrDropbox(locator.getPath())) {
                return null;   // Can't infer google & dropbox paths
            } else {
                String indexExtension =
                        (locator.getURLPath().toLowerCase().endsWith(".gz") || locator.getPath().toLowerCase().endsWith(".bgz")) ? ".tbi" : Tribble.STANDARD_INDEX_EXTENSION;
                return appendToPath(locator, indexExtension);
            }
        }
    }

    public void setAttribute(String key, Object value) {
        this.attributes.put(key, value);
    }

    public Object getAttribute(String key) {
        return attributes.get(key);
    }

    public void setIndexed(boolean indexed) {
        this.indexed = indexed;
    }

    public boolean isIndexed() {
        return indexed;
    }

    private static boolean isCloudOrDropbox(String path) {
        try {
            if (GoogleUtils.isGoogleDrive(path)) {
                return true;
            }
            if (path.startsWith("http://") || path.startsWith("https://")) {
                String host = new URL(path).getHost();
                if (host.equals("www.dropbox.com") || host.equals("dl.dropboxusercontent.com")) {
                    return true;
                }
            }
            return false;
        } catch (MalformedURLException e) {
            return false;
        }
    }

    /**
     * FOR LOAD FROM SERVER
     */
    public enum AttributeType {

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
        MAPPING("mapping"),
        COLOR("color"),
        INDEX("index");

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
