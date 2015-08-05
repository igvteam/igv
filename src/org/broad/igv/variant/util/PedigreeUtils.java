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
