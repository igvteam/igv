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

package org.broad.igv.variant.util;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * @author jrobinso
 * @date Apr 25, 2011
 */
public class PedigreeUtils {


    /**
     * Parses a pedigree file containing trios of child, father, mother.
     *
     * @param path
     * @throws IOException
     */

    public static void parseTrioFile(String path) throws IOException {

        BufferedReader reader = null;


        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine;
            int familyNumber = 1;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length == 3) {
                    String family = String.valueOf(familyNumber);
                    for (String member : tokens) {
                        AttributeManager.getInstance().addAttribute(member, "Family", family);
                    }
                    String child = tokens[0];
                    String father = tokens[1];
                    String mother = tokens[2];
                    AttributeManager.getInstance().addAttribute(child, "Member", "child");
                    AttributeManager.getInstance().addAttribute(child, "Relation", "child");
                    AttributeManager.getInstance().addAttribute(father, "Member", "father");
                    AttributeManager.getInstance().addAttribute(father, "Relation", "parent");
                    AttributeManager.getInstance().addAttribute(mother, "Member", "mother");
                    AttributeManager.getInstance().addAttribute(mother, "Relation", "parent");
                    familyNumber++;
                }
            }

            IGV.getInstance().notifyGroupEvent();
        } finally {
            if (reader != null) reader.close();
        }

    }
}
