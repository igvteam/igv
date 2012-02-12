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

package org.broad.igv.track;

import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Used to wrap multiple files as a single logical unit for tracks that break data by chromosome.
 */

public class MultiFileWrapper {

    String name;
    Map<String, ResourceLocator> locators;

    public MultiFileWrapper(String name, Map<String, ResourceLocator> locators) {
        this.name = name;
        this.locators = locators;
    }

    public String getName() {
        return name;
    }

    public Map<String, ResourceLocator> getLocators() {
        return locators;
    }

    public static MultiFileWrapper parse(ResourceLocator rl) {

        AsciiLineReader reader = null;

        try {
            reader = ParsingUtils.openAsciiReader(rl);

            String name = reader.readLine();
            String serverURL = reader.readLine();

            Map<String, ResourceLocator> locators = new LinkedHashMap();

            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];
                String file = tokens[1];
                ResourceLocator loc = new ResourceLocator(serverURL, file);
                locators.put(chr, loc);
            }

            return new MultiFileWrapper(name, locators);

        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (reader != null) {
                reader.close();

            }
        }

    }

}
