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

package org.broad.igv.feature.genome;

import java.io.IOException;
import java.io.InputStream;

public class GenomeResourceDescriptor extends GenomeDescriptor {

    public GenomeResourceDescriptor(String name,
                                    int version,
                                    String id,
                                    String cytoBandFileName,
                                    String geneFileName,
                                    String geneTrackName,
                                    String chrAliasFileName,
                                    String sequenceLocation,
                                    boolean chromosomesAreOrdered) {
        super(name, version, false, id, cytoBandFileName, geneFileName, chrAliasFileName,
                geneTrackName, sequenceLocation, chromosomesAreOrdered, false);
    }

    public InputStream getCytoBandStream() throws IOException {
        return GenomeResourceDescriptor.class.getResourceAsStream("/resources/hg18_cytoBand.txt");
        //return GenomeResourceDescriptor.class.getResourceAsStream(cytoBandFileName);
    }

    public InputStream getGeneStream() throws IOException {
        if (geneFileName != null) {
            return GenomeResourceDescriptor.class.getResourceAsStream(geneFileName);
        } else {
            return null;
        }
    }

    @Override
    public InputStream getChrAliasStream() throws IOException {
        if (chrAliasFileName != null) {
             return GenomeResourceDescriptor.class.getResourceAsStream(chrAliasFileName);
         } else {
             return null;
         }
    }

    public boolean isUserDefined() {
        return false;
    }

}
