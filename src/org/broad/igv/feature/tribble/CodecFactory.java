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

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.peaks.PeakCodec;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.util.BlockCompressedInputStream;
import org.broad.tribble.util.SeekableStreamFactory;
import org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFCodec;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * A factory class for Tribble codecs.  implements a single, static, public method to return a codec given a
 * path to a feature file  (bed, gff, vcf, etc).
 *
 */
public class CodecFactory {

    private static Logger log = Logger.getLogger(CodecFactory.class);

    /**
     * Return a tribble codec to decode the supplied file.
     *
     * @param path the path (file or URL) to the feature rile.
     */
    public static FeatureCodec getCodec(String path) {

        String fn = path.toLowerCase();
        if (fn.endsWith(".gz")) {
            int l = fn.length() - 3;
            fn = fn.substring(0, l);
        }

        if (fn.endsWith(".vcf4")) {
            return new VCFCodec();
        } else if (fn.endsWith(".vcf")) {
            return getVCFCodec(path);
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

        } else {
            throw new DataLoadException("Unknown file type", path);
        }

    }


    /**
     * Return the appropriate VCFCodec based on the version tag.
     * 
     * e.g.  ##fileformat=VCFv4.1
     *
     * @param path Path to the VCF file.  Can be a file path, or URL
     * @return
     */
    private static FeatureCodec getVCFCodec(String path) {

        BufferedReader reader = null;

        try {
            // If the file ends with ".gz" assume it is a tabix indexed file
            if (path.toLowerCase().endsWith(".gz")) {
                reader = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(
                        SeekableStreamFactory.getStreamFor(path))));
            } else {
                reader = ParsingUtils.openBufferedReader(path);
            }
            // Look for fileformat directive.  This should be the first line, but just in case check the first 20
            int lineCount = 0;
            String formatLine;
            while ((formatLine = reader.readLine()) != null && lineCount < 20) {
                if (formatLine.toLowerCase().startsWith("##fileformat")) {
                    String[] tmp = formatLine.split("=");
                    if (tmp.length > 1) {
                        String version = tmp[1].toLowerCase();
                        if (version.startsWith("vcfv3")) {
                            return new VCF3Codec();
                        } else {
                            return new VCFCodec();
                        }
                    }
                }
                lineCount++;
            }

        }
        catch (IOException e) {
            log.error("Error checking VCF Version");

        }
        finally {
            if (reader != null) try {
                reader.close();
            } catch (IOException e) {

            }
        }
        // Should never get here, but as a last resort assume this is a VCF 4.x file.
        return new VCFCodec();
    }
}
