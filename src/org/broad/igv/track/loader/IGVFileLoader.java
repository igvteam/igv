/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

package org.broad.igv.track.loader;

import org.broad.igv.data.DatasetDataSource;
import org.broad.igv.data.IGVDataset;
import org.broad.igv.data.IGVDatasetParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.ScatterplotRenderer;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;

import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 * @date Dec 16, 2010
 */
public class IGVFileLoader extends AbstractLoader implements TrackLoader {


    public boolean canLoad(ResourceLocator locator) {

        String typeString = getTypeString(locator);

        if(typeString.endsWith(".cn") || typeString.endsWith("xcn") || typeString.endsWith("snp") ||
                typeString.endsWith(".igv") || typeString.endsWith("loh")) {
            return true;
        }
        try {
            return IGVDatasetParser.parsableMAGE_TAB(locator);
        }
        catch(Exception e) {
            // Nothing to do here, apparently not a mage tab file
        }
        return false;
    }

    /**
     * Load the resource and return the newly created tracks
     * @param locator
     * @return
     */
    public List<Track> load(ResourceLocator locator) {

        //if (locator.isLocal()) {
        //    if (!checkSize(locator.getPath())) {
        //        return;
        //    }
        //}

        List<Track> newTracks = new ArrayList();

        String dsName = locator.getTrackName();
        String currentGenomeId = GenomeManager.getInstance().getGenomeId();

        IGVDataset ds = new IGVDataset(currentGenomeId, locator);
        ds.setName(dsName);

        TrackProperties trackProperties = ds.getTrackProperties();
        String path = locator.getPath();
        TrackType type = ds.getType();

        for (String trackName : ds.getTrackNames()) {

            Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
            DatasetDataSource dataSource = new DatasetDataSource(currentGenome, trackName, ds);
            String trackId = path + "_" + trackName;
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            // track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());
            track.setTrackProperties(trackProperties);

            if (type == TrackType.ALLELE_FREQUENCY) {
                track.setRendererClass(ScatterplotRenderer.class);
                track.setHeight(40);
            }
            newTracks.add(track);
        }

        return newTracks;
    }
}
