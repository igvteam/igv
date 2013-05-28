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

import org.broad.igv.cli_plugin.LineFeatureDecoder;

/**
 * Wrapper class to make it easy to just extract the first sample from
 * an fpkm_tracking file.
 * @author jacob
 * @date 2013-May-28
 */
public class FPKMTrackingSampleCodec extends CufflinksCodec<FPKMSampleValue>  implements LineFeatureDecoder<FPKMSampleValue> {

    private FPKMTrackingCodec trackingCodec;

    public FPKMTrackingSampleCodec() {
        super(FPKMSampleValue.class, "Plugin");
        this.trackingCodec = new FPKMTrackingCodec(this.path);
    }

    @Override
    protected Object readHeader(String[] tokens) {
        return trackingCodec.readHeader(tokens);
    }

    @Override
    public FPKMSampleValue decode(String line) {
        return trackingCodec.decode(line).getSampleValue(0);
    }
}
