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

package org.broad.igv.gs.dm;


import java.io.File;
import java.net.URL;
import java.util.List;

/**
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSDirectoryListing {

    private GSFileMetadata directory;

    private List<GSFileMetadata> contents;

    public GSDirectoryListing(String url, List<GSFileMetadata> contents) {

        // Parse URL to get components.
        //URL u = new URL(url);
        String path = null;
        try {
            path = (new URL(url)).getPath();
        } catch (Exception e) {

        }
        String name = (new File(path)).getName();
        String format = "";
        String size = "";
        boolean isDirectory = true;
        directory = new GSFileMetadata(name, path, url, format, size, isDirectory);
        this.contents = contents;
    }

    public GSFileMetadata getDirectory() {
        return directory;
    }

    public List<GSFileMetadata> getContents() {
        return contents;
    }
}
