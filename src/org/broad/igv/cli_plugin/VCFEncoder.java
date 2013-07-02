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

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFVariant;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFHeader;

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
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<VCFVariant> features) throws IOException {

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
