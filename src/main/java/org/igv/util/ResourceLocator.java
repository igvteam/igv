package org.igv.util;

import htsjdk.tribble.Tribble;
import org.igv.feature.genome.load.TrackConfig;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.TrackProperties;
import org.json.JSONObject;

import java.awt.*;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.*;

import static org.igv.feature.tribble.CodecFactory.ucscSNP;

//import java.awt.*;

/**
 * Represents a data file or other resource, which might be local file or remote resource.
 *
 * @author jrobinso
 */
public class ResourceLocator {

    private static final Logger log = LogManager.getLogger(ResourceLocator.class);

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

    String labelField;

    /**
     * Descriptive text
     */
    String description;

    /**
     * The type of resource (generally this refers to the file format)
     */
    public String format;

    /**
     * A UCSC style track line.  Overrides value in file, if any.
     */
    String trackLine;  //

    TrackProperties trackProperties;

    /**
     * Color for features or data.  Somewhat redundant with trackLine.
     */
    Color color;

    String sampleId;


    /**
     * Track metadata.  Primarily for generating description and popup text.
     */
    private Map<String, String> metadata;
    private boolean indexed;


    /**
     * True if this is an htsget resource
     */
    private boolean htsget;
    private boolean dataURL;
    private Integer visibilityWindow;
    private String trixURL;
    private String panelName;
    private String autoscaleGroup;
    private String[] filterTypes;

    public static List<ResourceLocator> getLocators(Collection<File> files) {

        List<ResourceLocator> locators = new ArrayList<>();

        Set<String> indexExtensions = new HashSet<>(Arrays.asList("bai", "crai", "sai", "tbi", "tbx"));
        Set<File> indexes = new HashSet<>();
        Map<String, File> indexMap = new HashMap<>();
        for (File f : files) {
            String fn = f.getName();
            int idx = fn.lastIndexOf('.');
            if (idx > 0) {
                String ext = fn.substring(idx + 1);
                if (indexExtensions.contains(ext)) {
                    String base = fn.substring(0, idx);
                    if (ext.equals(".bai") && !base.endsWith(".bam")) {
                        base += ".bam";   // Picard convention
                    } else if (ext.equals(".crai") && !base.endsWith(".cram")) {
                        base += ".cram";  // Possible Picard convention
                    }
                    indexes.add(f);
                    indexMap.put(base, f);
                }
            }
        }

        for (File f : files) {
            if (indexes.contains(f)) continue;
            ResourceLocator locator = new ResourceLocator(f.getAbsolutePath());
            File indexFile = indexMap.get(f.getName());
            if (indexFile != null) {
                locator.setIndexPath(indexFile.getAbsolutePath());
            }
            locators.add(locator);
        }

        return locators;
    }


    /**
     * Create a ResourceLocator from a TrackConfig object.  This is used to create
     * a ResourceLocator for a track hub track.
     *
     * @param trackConfig
     * @return
     */

    public static ResourceLocator fromTrackConfig(TrackConfig trackConfig) {
        String trackPath = trackConfig.url;
        ResourceLocator res = new ResourceLocator(trackPath);
        res.setName(trackConfig.name);
        res.setIndexPath(trackConfig.indexURL);
        res.setFormat(trackConfig.format);
        res.setVisibilityWindow(trackConfig.visibilityWindow);
        res.setPanelName(trackConfig.panelName);
        res.setTrixURL(trackConfig.trixURL);
        res.setFeatureInfoURL(trackConfig.infoURL);
        res.setIndexed(trackConfig.indexed != null ? trackConfig.indexed : false);
        res.setLabelField(trackConfig.labelField);
        res.setDescription(trackConfig.description);
        res.setAutoscaleGroup(trackConfig.autoscaleGroup);
        res.setTrackProperties(new TrackProperties(trackConfig));
        res.setFilterTypes(trackConfig.filterTypes);
        return res;

    }

    private void setFilterTypes(String[] filterTypes) {
        this.filterTypes = filterTypes;
    }

    public String[] getFilterTypes() {
        return filterTypes;
    }

    /**
     * Constructor for local files and URLs
     *
     * @param path
     */
    public ResourceLocator(String path) {

        this.setPath(path);
        if (path != null && path.startsWith("https://") && GoogleUtils.isGoogleDrive(path)) {
            this.resolveGoogleDrive(path);
        } else if (path != null && path.startsWith("htsget://")) {
            this.htsget = true;
        }

    }

    private void resolveGoogleDrive(String path) {

        JSONObject fileInfo = GoogleUtils.getDriveFileInfo(path);
        this.name = fileInfo.optString("name", null);
        this.format = deriveFormat(this.name);
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

    public boolean isDataURL() {
        return dataURL;
    }

    /**
     * Determines if the resource actually exists.
     *
     * @return true if resource was found.
     */
    public boolean exists() {
        return ParsingUtils.fileExists(path);
    }


    public void setFormat(String formatOrExt) {
        this.format = formatOrExt == null ? null : formatOrExt.startsWith(".") ? formatOrExt.substring(1) : formatOrExt;
    }

    /**
     * Return a string suitable for determining file type based on extension
     * May or may not be a full, readable path. txt and gz extensions are stripped
     *
     * @return
     */
    public String getFormat() {
        if (format == null) {
            format = deriveFormat(this.path);
        }
        return format;
    }


    public Integer getVisibilityWindow() {
        return visibilityWindow;
    }

    public void setVisibilityWindow(Integer visibilityWindow) {
        this.visibilityWindow = visibilityWindow;
    }

    public static String deriveFormat(String pathOrName) {

        String filename;
        if (FileUtils.isRemote(pathOrName)) {
            try {
                URL url = HttpUtils.createURL(pathOrName);
                filename = url.getPath().toLowerCase();
                String query = url.getQuery();
                if (query != null) {
                    Map<String, String> queryMap = URLUtils.parseQueryString(query);
                    // If type is set explicitly use it
                    if (queryMap.containsKey("dataformat")) {
                        return queryMap.get("dataformat");
                    } else if (queryMap.containsKey("file")) {
                        filename = queryMap.get("file");
                    }
                }

            } catch (MalformedURLException e) {
                log.error("Error interpreting url: " + pathOrName, e);
                filename = pathOrName;
            }
        } else {
            filename = pathOrName.toLowerCase();
        }

        // Strip .txt, .gz, and .xls extensions.  (So  foo.cn.gz => a cn file)
        if ((filename.endsWith(".txt") || filename.endsWith(".xls") || filename.endsWith(".gz") ||
                filename.endsWith(".bgz"))) {
            filename = filename.substring(0, filename.lastIndexOf(".")).trim();
        }

        // Some special cases
        if (filename.endsWith("fpkm_tracking")) {
            return "fpkm_tracking";
        } else if (filename.endsWith("gene_exp.diff")) {
            return "gene_exp.diff";
        } else if (filename.endsWith("cds_exp.diff")) {
            return "cds_exp.diff";
        } else if (filename.endsWith("genepredext")) {
            return "genepredext";
        } else if (filename.endsWith(".maf.annotated")) {
            // TCGA extension
            return "mut";
        } else if (filename.endsWith("junctions.bed")) {
            return "junctions";
        } else if (filename.endsWith("bam.list")) {
            return "bam.list";
        } else if (filename.endsWith("sam.list")) {
            return "sam.list";
        } else if (filename.endsWith("vcf.list")) {
            return "vcf.list";
        } else if (filename.endsWith("seg.zip")) {
            return "seg.zip";
        } else if (filename.endsWith(".ewig.tdf") || (filename.endsWith(".ewig.ibf"))) {
            return "ewig.tdf";
        } else if (filename.endsWith(".maf.dict")) {
            return "maf.dict";
        } else if (filename.endsWith("_clusters")) {
            return "bedpe";
        } else {
            // Default - derive format from the extension.  If not a common format check special cases
            String format = filename.substring(filename.lastIndexOf('.') + 1).toLowerCase();
            if (!knownFormats.contains(format)) {
                if (filename.contains("refflat")) {
                    return "reflat";
                } else if (filename.contains("genepred") || filename.contains("ensgene") ||
                        filename.contains("refgene") || filename.contains("ncbirefseq")) {
                    return "refgene";
                } else if (filename.contains("ucscgene")) {
                    return "ucscgene";
                } else if (filename.matches(ucscSNP)) {
                    return "snp";
                }
            }
            return format;
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
        if (path.startsWith("http://") || path.startsWith("https://") || path.startsWith("gs://")) {
            int idxQuestion = path.indexOf('?');
            String actualPath = idxQuestion < 0 ? path : path.substring(0, idxQuestion);
            int idxSlash = actualPath.lastIndexOf('/');
            return idxSlash < 0 ?
                    actualPath :
                    actualPath.substring(idxSlash + 1);

        } else {
            return (new File(path)).getName();
        }
    }


    public String getDBUrl() {
        return dbURL;
    }

    public boolean isLocal() {
        return dbURL == null && !dataURL && !FileUtils.isRemote(path);
    }

    public void setTrackInfoURL(String trackInforURL) {
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
        if (name != null) {
            return name;
        } else {
            return this.getFileName();
        }
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public void setTrackProperties(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
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

    public String getLabelField() {
        return labelField;
    }

    public void setLabelField(String labelField) {
        this.labelField = labelField;
    }

    public void setFeatureInfoURL(String featureInfoURL) {
        this.featureInfoURL = featureInfoURL;
    }

    public void setPath(String path) {

        if (path != null && path.startsWith("file://")) {
            this.path = path.substring("file://".length());
        } else if (path != null && path.startsWith("s3://")) {
            this.path = path;
            String s3UrlIndexPath = detectIndexPath(path);
            this.setIndexPath(s3UrlIndexPath);
        } else {
            this.dataURL = ParsingUtils.isDataURL(path);
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
     * Add the {@code extension} to the path in locator, preserving
     * query string elements if present
     *
     * @param locator
     * @param extension
     * @return
     */
    public static String appendToPath(ResourceLocator locator, String extension) {
        String extendedPath = locator.getURLPath() + extension;
        String qs = locator.getURLQueryString();
        if (qs != null && qs.length() > 0) {
            extendedPath += "?" + qs;
        }
        return extendedPath;
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

    public void setIndexed(boolean indexed) {
        this.indexed = indexed;
    }

    public boolean isIndexed() {
        return indexed;
    }

    public boolean isHtsget() {
        return htsget;
    }

    public void setHtsget(boolean htsget) {
        this.htsget = htsget;
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

    public void setTrixURL(String trixURL) {
        this.trixURL = trixURL;
    }

    public String getTrixURL() {
        return trixURL;
    }

    public Map<String, String> getMetadata() {
        return metadata;
    }

    public void setMetadata(Map<String, String> metadata) {
        this.metadata = metadata;
    }

    public String getPanelName() {
        return panelName;
    }

    public void setPanelName(String panelName) {
        this.panelName = panelName;
    }

    public void setDataURL(boolean dataURL) {
        this.dataURL = dataURL;
    }

    public boolean isDataUrl() {
        return dataURL;
    }

    public void setAutoscaleGroup(String autoscaleGroup) {
        this.autoscaleGroup = autoscaleGroup;
    }

    public String getAutoscaleGroup() {
        return autoscaleGroup;
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
        LABEL_FIELD("labelField"),
        RESOURCE_TYPE("resourceType"),
        TRACK_LINE("trackLine"),
        COVERAGE("coverage"),
        MAPPING("mapping"),
        COLOR("color"),
        INDEX("index"),
        HTSGET("htsget");

        private final String name;

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

    static Set<String> knownFormats = new HashSet<>(Arrays.asList("gff", "bed", "gtf", "gff3",
            "seg", "bb", "bigbed", "bigwig", "bam", "cram", "vcf", "bedmethyl"));
}
