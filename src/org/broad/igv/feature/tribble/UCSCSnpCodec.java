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

import org.broad.igv.Globals;
import org.broad.igv.feature.UCSCSnpFeature;
import org.broad.igv.feature.genome.Genome;

/**
 * Created by jrobinso on 5/26/15.
 */
public class UCSCSnpCodec extends UCSCCodec<UCSCSnpFeature> {

    private Genome genome;

    public UCSCSnpCodec(Genome genome) {
        super(UCSCSnpFeature.class);
        this.genome = genome;
    }

    @Override
    public UCSCSnpFeature decode(String s) {

        String[] tokens = Globals.tabPattern.split(s);

        if (tokens.length < 25) return null;

        String chr = tokens[1];
        if (genome != null) {
            chr = genome.getChromosomeAlias(chr);
        }
        int start = Integer.parseInt(tokens[2]);
        int end = Integer.parseInt(tokens[3]);
        return new UCSCSnpFeature(chr, start, end, tokens);

    }

    @Override
    public boolean canDecode(String path) {
        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.toLowerCase().endsWith(".snp");
    }
}
