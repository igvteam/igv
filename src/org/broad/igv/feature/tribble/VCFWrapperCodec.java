/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.feature.tribble;

import org.broad.igv.variant.vcf.VCFVariant;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * @author Jim Robinson
 * @date Aug 1, 2011
 */
public class VCFWrapperCodec implements FeatureCodec {

    FeatureCodec wrappedCodec;

    public VCFWrapperCodec(FeatureCodec wrappedCodec) {
       this.wrappedCodec = wrappedCodec;
    }

    public Feature decodeLoc(String line) {
        return wrappedCodec.decodeLoc(line);
    }

    public Feature decode(String line) {
        VariantContext vc = (VariantContext) wrappedCodec.decode(line);
        return vc == null ? null : new VCFVariant(vc);
    }

    public Class getFeatureType() {
        return VCFVariant.class;
    }

    public Object readHeader(LineReader reader) {
        return wrappedCodec.readHeader(reader);
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        return path.endsWith(".vcf") || path.endsWith(".vcf4");
    }
}
