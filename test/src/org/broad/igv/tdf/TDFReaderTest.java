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

package org.broad.igv.tdf;

import org.broad.igv.util.ResourceLocator;
import org.junit.Test;

import java.util.Set;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 *         Date: 7/19/12
 *         Time: 11:29 AM
 */
public class TDFReaderTest {

    @Test
    public void testReader() throws Exception {
        String url = "http://www.broadinstitute.org/igvdata/encode/hg18/broadHistone/SignalK562H3k4me3.tdf";

        TDFReader reader =  new TDFReader(new ResourceLocator(url));

        int version = reader.getVersion();
        assertEquals(3, version);
        assertTrue(reader.compressed);

        String [] trackNames = reader.getTrackNames();
        int nTracks = trackNames.length;
        assertEquals(1, nTracks);

        String trackName = trackNames[0];
        assertEquals(trackName, "SignalK562H3k4me3.wig.gz");

        Set<String> chrNames = reader.getChromosomeNames();
        int nChromosomes = chrNames.size();
        assertEquals(24, nChromosomes);

        String datasetName =  "/chr1/z0/mean";
        TDFDataset dataset = reader.getDataset(datasetName);
        assertEquals(datasetName, dataset.getName());

        TDFTile tile = reader.readTile(dataset, 0);
        assertNotNull(tile);

    }
}
