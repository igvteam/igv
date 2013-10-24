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

package org.broad.igv.feature.tribble.reader;

import org.broad.tribble.readers.*;

import java.io.IOException;


import java.io.IOException;

/**
 * @author Jim Robinson
 * @date 2/11/12
 */
public class TabixIteratorLineReader implements LineReader {

    TabixReader.Iterator iterator;


    public TabixIteratorLineReader(TabixReader.Iterator iterator) {
        this.iterator = iterator;
    }

    public String readLine() {
        try {
            return iterator != null ? iterator.next() : null;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void close() {
        // Ignore -
    }
}
