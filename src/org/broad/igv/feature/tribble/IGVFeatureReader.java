/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.tribble;




import htsjdk.tribble.Feature;

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
