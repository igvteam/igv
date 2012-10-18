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

import java.io.IOException;
import java.io.OutputStream;
import java.util.Iterator;
import java.util.Map;

/**
 * Interface for creating a string from a feature.
 * User: jacob
 * Date: 2012-Aug-02
 *
 * @see PluginCodec
 * @see FeatureDecoder
 */
public interface FeatureEncoder<T extends Feature> extends FeatureIO {


    /**
     * Write all of the {@code features} to {@code outputStream}
     *
     * @param outputStream
     * @param features
     * @return A map containing any attributes deemed necessary. This map will be provided
     *         to the {@link FeatureDecoder}. It may be null
     */
    Map<String, Object> encodeAll(OutputStream outputStream, Iterator<T> features) throws IOException;
}
