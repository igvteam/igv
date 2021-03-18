package org.broad.igv.sam.reader;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.cram.IGVReferenceSource;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.URLUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

/**
 * Simple pool for reusing SamReader instances.  The SamReader query object is not thread safe, so if triggering
 * multiple queries in parallel, which can easily occur in IGV,  a new SamReader instance is needed for each query.
 * When the query is completed "freeReader(reader)" is called which makes the reader available for future queries.
 * Its assumed that SamReader is lightweight, no attempt is made to limit the pool size or dispose of old instances.
 * In practice the pool will rarely grow to more than a few readers.
 */
public class SamReaderPool {

    private static Logger log = Logger.getLogger(SamReaderPool.class);

    private ResourceLocator locator;
    private boolean requireIndex;
    private List<SamReader> availableReaders;

    public SamReaderPool(ResourceLocator locator, boolean requireIndex) {
        this.locator = locator;
        this.requireIndex = requireIndex;
        availableReaders = Collections.synchronizedList(new ArrayList<>());
    }

    public SamReader getReader() throws IOException {
        System.out.println("Pool size = " + availableReaders.size());
        if (availableReaders.size() > 0) {
            return availableReaders.remove(0);
        } else {
            return createReader();
        }
    }

    public void freeReader(SamReader reader) {
        availableReaders.add(reader);
    }

    public void close() throws IOException {
        for(SamReader reader : availableReaders) {
            reader.close();
        }
        availableReaders.clear();
    }

    private SamReader createReader() throws IOException {

        boolean isLocal = locator.isLocal();
        final SamReaderFactory factory = SamReaderFactory.makeDefault().
                referenceSource(new IGVReferenceSource()).
                validationStringency(ValidationStringency.SILENT);
        SamInputResource resource;

        if (isLocal) {
            resource = SamInputResource.of(new File(locator.getPath()));
        } else {
            URL url = HttpUtils.createURL(locator.getPath());
            if (requireIndex) {
                resource = SamInputResource.of(IGVSeekableStreamFactory.getInstance().getStreamFor(url));
            } else {
                resource = SamInputResource.of(new BufferedInputStream(HttpUtils.getInstance().openConnectionStream(url)));
            }
        }

        if (requireIndex) {

            String indexPath = getExplicitIndexPath(locator);
            if (indexPath == null || indexPath.length() == 0) {
                indexPath = getIndexPath(locator.getPath());
            }
            if (isLocal) {
                File indexFile = new File(indexPath);
                resource = resource.index(indexFile);
            } else {
                SeekableStream indexStream = IGVSeekableStreamFactory.getInstance().getStreamFor(HttpUtils.createURL(indexPath));
                resource = resource.index(indexStream);
            }
        }

        return factory.open(resource);
    }

    /**
     * Fetch an explicitly set index path, either via the ResourceLocator or as a parameter in a URL
     *
     * @param locator
     * @return the index path, or null if no index path is set
     */
    private String getExplicitIndexPath(ResourceLocator locator) {

        String p = locator.getPath().toLowerCase();
        String idx = locator.getIndexPath();

        if (idx == null && (p.startsWith("http://") || p.startsWith("https://"))) {
            try {
                URL url = HttpUtils.createURL(locator.getPath());
                String queryString = url.getQuery();
                if (queryString != null) {
                    Map<String, String> parameters = URLUtils.parseQueryString(queryString);
                    if (parameters.containsKey("index")) {
                        idx = parameters.get("index");

                    }
                }
            } catch (MalformedURLException e) {
                log.error("Error parsing url: " + locator.getPath());
            }
        }
        return idx;
    }

    /**
     * Try to guess the index path.
     *
     * @param pathOrURL
     * @return
     * @throws IOException
     */
    private String getIndexPath(String pathOrURL) throws IOException {

        List<String> pathsTried = new ArrayList<String>();

        String indexPath;

        if (URLUtils.isURL(pathOrURL)) {

            String path = URLUtils.getPath(pathOrURL);

            if (path.endsWith(".bam")) {
                // Try .bam.bai
                indexPath = URLUtils.addExtension(pathOrURL, ".bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .bai
                indexPath = URLUtils.replaceExtension(pathOrURL, ".bam", ".bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .bam.csi
                indexPath = URLUtils.addExtension(pathOrURL, ".csi");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .csi
                indexPath = URLUtils.replaceExtension(pathOrURL, ".bam", ".csi");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }
            }

            //  cram
            if (path.endsWith(".cram")) {
                indexPath = URLUtils.addExtension(pathOrURL, ".crai");
                if (FileUtils.resourceExists(indexPath)) {
                    pathsTried.add(indexPath);
                    return indexPath;
                } else {
                    indexPath = pathOrURL.substring(0, pathOrURL.length() - 5) + ".crai";
                    if (FileUtils.resourceExists(indexPath)) {
                        return indexPath;
                    }
                }
            }


        } else {
            // Local file

            indexPath = pathOrURL + ".bai";

            if (FileUtils.resourceExists(indexPath)) {
                return indexPath;
            }

            if (indexPath.contains(".bam.bai")) {
                indexPath = indexPath.replaceFirst(".bam.bai", ".bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            } else {
                indexPath = indexPath.replaceFirst(".bai", ".bam.bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            }


            // Try .bam.csi
            indexPath = pathOrURL + ".csi";
            pathsTried.add(indexPath);
            if (FileUtils.resourceExists(indexPath)) {
                return indexPath;
            }

            // Try .csi
            if (pathOrURL.endsWith(".bam")) {
                indexPath = pathOrURL.substring(0, pathOrURL.length() - 4) + ".csi";
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            }

            if (pathOrURL.endsWith(".cram")) {
                indexPath = pathOrURL + ".crai";
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                } else {
                    indexPath = pathOrURL.substring(0, pathOrURL.length() - 5) + ".crai";
                    if (FileUtils.resourceExists(indexPath)) {
                        return indexPath;
                    }
                }
            }


        }

        String defaultValue = pathOrURL + (pathOrURL.endsWith(".cram") ? ".crai" : ".bai");
        indexPath = MessageUtils.showInputDialog(
                "Index is required, but no index found.  Please enter path to index file:",
                defaultValue);
        if (indexPath != null && FileUtils.resourceExists(indexPath)) {
            return indexPath;
        }


        String msg = "Index file not found.  Tried ";
        for (String p : pathsTried) {
            msg += "<br>" + p;
        }
        throw new DataLoadException(msg, indexPath);

    }

}
