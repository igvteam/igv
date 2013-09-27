/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.api;

import org.broad.igv.track.Track;

import java.io.IOException;
import java.util.List;

/**
 * Handler for loading data from a file or URL.
 *
 * @author jacob
 * @since 2013-Sep-27
 * @api
 */
public interface LoadHandler {

    /**
     * Called when a file at {@code path} is to be loaded.
     * {@code path} may be a file or URL. The generated
     * tracks should be added to {@code newTracks}
     * @param path
     * @param newTracks
     * @throws IOException
     */
    public void load(String path, List<Track> newTracks)  throws IOException;
}
