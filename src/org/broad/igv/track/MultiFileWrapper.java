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

package org.broad.igv.track;

import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

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
                ResourceLocator loc = new ResourceLocator(file);
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
