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

package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.Range;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.AsciiFeatureCodec;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Parses various cufflinks (and cuffdiff) output files as described here:
 * <p/>
 * http://cufflinks.cbcb.umd.edu/manual.html
 *
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:30 PM
 */
public class CufflinksParser {

    public static List<? extends Range> parse(String path) throws IOException {

        final String s = path.toLowerCase();
        if (s.endsWith("fpkm_tracking")) {
            AsciiFeatureCodec<FPKMValue> codec = new FPKMTrackingCodec(path);
            return parse(codec, path);
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            AsciiFeatureCodec<ExpDiffValue> codec = new ExpDiffCodec(path);
            return parse(codec, path);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }

    }

    public static <T extends Range> List<T> parse(AsciiFeatureCodec<T> codec, String path) throws IOException {
        List<T> values = new ArrayList<T>();
        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(path, codec, false);
        Iterator<T> iter = reader.iterator();
        while(iter.hasNext()){
            values.add(iter.next());
        }
        return values;
    }

}
