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

package org.broad.igv.cli_plugin;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFHeader;

import java.io.IOException;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @date 2013-Jul-02
 */
public class VCFEncoder implements FeatureEncoder<VCFVariant> {

    private VCFHeader header;

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        assert argument.getType() == Argument.InputType.VARIANT_TRACK;
        VariantTrack track = (VariantTrack) argumentMap.get(argument);
        this.header = (VCFHeader) track.getHeader();
    }

    @Override
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<? extends VCFVariant> features) throws IOException {

        SAMSequenceDictionary seqDict = new SAMSequenceDictionary();
        EnumSet<Options> options = VariantContextWriterFactory.DEFAULT_OPTIONS;
        options.add(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        options.remove(Options.INDEX_ON_THE_FLY);
        VariantContextWriter writer = VariantContextWriterFactory.create(outputStream, seqDict, options);

        writer.writeHeader(this.header);

        while(features.hasNext()){
            writer.add(features.next().getVariantContext());
        }

        return null;
    }
}
