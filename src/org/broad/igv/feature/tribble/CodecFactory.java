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

import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.peaks.PeakCodec;
import org.broad.tribble.FeatureCodec;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 28, 2010
 * Time: 11:40:43 AM
 * To change this template use File | Settings | File Templates.
 */
public class CodecFactory {
    public static FeatureCodec getCodec(String file) {

        String fn = file.toLowerCase();
        if (fn.endsWith(".gz")) {
            int l = fn.length() - 3;
            fn = fn.substring(0, l);
        }
        if (fn.endsWith(".vcf") || fn.endsWith(".vcf4") ) {
            return new org.broad.tribble.vcf.VCFCodec();
        } else if (fn.endsWith(".bed")) {
            return new BEDCodec();
        } else if (fn.endsWith(".repmask")) {
            return new REPMaskCodec();
        } else if (fn.endsWith(".gff3")) {
            return new GFFCodec(GFFCodec.Version.GFF3);
        } else if (fn.endsWith(".gff")) {
            return new GFFCodec();
        } else if (fn.endsWith(".sam")) {
            return new SAMCodec();
        } else if (fn.endsWith(".psl") || fn.endsWith(".pslx")) {
            return new PSLCodec();

        } else if (fn.endsWith(".peak")) {
            return new PeakCodec();

        }
        else {
            throw new DataLoadException("Unknown file type", file);
        }

    }
}
