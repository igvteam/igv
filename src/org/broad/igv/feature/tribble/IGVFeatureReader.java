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

package org.broad.igv.feature.tribble;




import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Interface to replace the Tribble equivalent.  Returns plain "Iterator" rather than "CloseableIterator".  A wrapper
 * implementation is provided that handles closing the tribble resource.
 *
 * NOTE:  IGV classes should use this interface, rather than the Tribble FeatureReader.
 *
 * @author jrobinso
 *         Date: 5/8/13
 *         Time: 10:55 AM
 */
public interface IGVFeatureReader {

    public Iterator<Feature> query(final String chr, final int start, final int end) throws IOException;

    public Iterator<Feature> iterator() throws IOException;

    public List<String> getSequenceNames();

    public Object getHeader();
}
