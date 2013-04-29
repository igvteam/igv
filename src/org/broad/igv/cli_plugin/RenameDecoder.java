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

import org.broad.igv.feature.BasicFeature;
import org.broad.tribble.AsciiFeatureCodec;

/**
 * The purpose of this class is to rename incoming features, if we don't like the
 * way the external tool names them. For instance, Cufflinks uses a numbering scheme
 * This doesn't work for on-the-fly calculation since each feature will always be
 * CUFF1.1 / 1.2 etc. We would like the naming to be consistent between calculations.
 *
 * @author jacob
 * @date 2013-Apr-29
 */
public class RenameDecoder<T extends BasicFeature> extends AsciiDecoder<T>{

    public RenameDecoder(AsciiFeatureCodec<T> featureCodec){
        super(new AsciiDecoder.DecoderWrapper<T>(featureCodec));
    }

    @Override
    public T decode(String line) {
        T bf = super.decode(line);
        bf.setName(createName(bf));
        return bf;
    }

    /**
     * Create a new name for this feature.
     * Intended to be overridden, default implementation uses start position.
     * @param bf
     * @return
     */
    protected String createName(BasicFeature bf) {
        return String.format("%d", bf.getStart());
    }
}
