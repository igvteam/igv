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

package org.broad.igv.ui.genome;

import java.util.HashMap;
import java.util.Map;

/**
 * A record representing a genome hosted on a remote server, or more specifically a row in the IGV hosted genome table.
 */
public class GenomeListItem {

    String id;
    String path;
    String displayableName;

    Map<String, String> attributes;

    public GenomeListItem(String displayableName, String path, String id, Map<String, String> attributes) {
        this.id = id;
        this.path = path;
        this.displayableName = displayableName;
        this.attributes = attributes;
    }

    /**
     * Alternate constructor for legacy genome list format.
     *
     * @param displayableName The name that can be shown to a user.
     * @param path            The location of the genome archive, can be a file path or URL
     * @param id              The id of the genome.
     */
    public GenomeListItem(String displayableName, String path, String id) {
        this.id = id;
        this.path = path;
        this.displayableName = displayableName;
        attributes = new HashMap<>();
        attributes.put("common name", displayableName);
        attributes.put("url", path);
        attributes.put("assembly", id);
    }

    public String getAttributeValue(String name) {
        return attributes.get(name);
    }

    public String getId() {
        return id;
    }

    public String getPath() {
        return path;
    }

    public String getDisplayableName() {
        return displayableName;
    }

    public String toString() {
        return displayableName;
    }

@Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GenomeListItem)) return false;
        GenomeListItem that = (GenomeListItem) o;
        return id != null && id.equals(that.id) &&
               path != null && path.equals(that.path);
    }

    @Override
    public int hashCode() {
        int result = id != null ? id.hashCode() : 0;
        result = 31 * result + (path != null ? path.hashCode() : 0);
        return result;
    }
}
