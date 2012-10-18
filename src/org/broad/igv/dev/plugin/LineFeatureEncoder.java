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

package org.broad.igv.dev.plugin;

import org.broad.tribble.Feature;

/**
 * User: jacob
 * Date: 2012-Oct-01
 */
public interface LineFeatureEncoder<T extends Feature> extends FeatureIO {

    /**
     * Create a single data line (excluding newline character)
     * from this Feature. Return null to skip encoding
     *
     * @param feature
     * @return
     */
    String encode(T feature);

    /**
     * @param line
     * @return The number of data columns contained in this line.
     *         Some decoders need this information.
     */
    int getNumCols(String line);

    /**
     * Get header for this encoder. This will
     * be written before any features.
     *
     * @return String header, or null for none
     */
    String getHeader();
}
