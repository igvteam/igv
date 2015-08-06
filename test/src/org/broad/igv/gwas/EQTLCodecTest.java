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

package org.broad.igv.gwas;

import org.broad.igv.feature.genome.Genome;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 4/26/13
 *         Time: 1:36 PM
 */
public class EQTLCodecTest {
    @Test
    public void testDecode() throws Exception {

        Genome genome = null;  // Genome not needed for this test

        EQTLCodec codec = new EQTLCodec(genome);

        String line = "rs4700772\t5\t180341845\tENSG00000168903.7\tBTNL3\t180415845\t129.858921455635\t1.61918925670937e-27";

        EQTLFeature f = codec.decode(line);

        assertEquals("rs4700772", f.getSnp());
        assertEquals("ENSG00000168903.7", f.getGeneId());
        assertEquals("BTNL3", f.getGeneName());
        assertEquals("5", f.getChr());
        assertEquals(180341844, f.getStart());
        assertEquals(180341845, f.getEnd());

    }
}
