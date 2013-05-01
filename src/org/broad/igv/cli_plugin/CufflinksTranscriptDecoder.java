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

package org.broad.igv.cli_plugin;

import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.track.GFFFeatureSource;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

/**
 * @author jacob
 * @date 2013-Apr-29
 */
public class CufflinksTranscriptDecoder extends RenameDecoder{
    public CufflinksTranscriptDecoder() {
        super(new GFFCodec(GenomeManager.getInstance().getCurrentGenome()));
    }

    @Override
    public Iterator decodeAll(InputStream is, boolean strictParsing) throws IOException {
        Iterator iter = super.decodeAll(is, strictParsing);
        GFFFeatureSource.GFFCombiner combiner = new GFFFeatureSource.GFFCombiner();
        combiner.addFeatures(iter);
        return combiner.combineFeatures().iterator();
    }
}
