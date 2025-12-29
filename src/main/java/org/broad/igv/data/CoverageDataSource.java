package org.broad.igv.data;

/**
 * @author Fabien Campagne
 *         Date: 6/10/11
 *         Time: 4:45 PM
 */
public interface CoverageDataSource extends DataSource {
    /**
     * The filename that contains the coverage data. Used to persist the state of coverage tracks.
     * @return a filename.
     */
    String getPath();

    /**
     * Tell the coverage source to normalize coverage by some appropriate normalization method.
     * @param normalize True if normalization should be performed, false otherwise.
     */
    void setNormalize(boolean normalize);

    /**
     * Whether the source is performing
     * some normalization operation
     */
    boolean getNormalize();
}
