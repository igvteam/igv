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

package org.broad.igv.feature.genome;

import org.broad.igv.feature.genome.fasta.FastaUtils;
import org.broad.igv.util.HttpUtils;

import java.io.File;
import java.io.IOException;
import java.util.Objects;
import java.util.zip.ZipException;

/**
 * A container for specific genome information which can be used to
 * manage loaded genomes.
 */
public class GenomeListItem {

    private String displayableName;
    private String path;
    private String id;
    private boolean hasDownloadedSequence = false;

    public static final GenomeListItem ITEM_MORE;

    static {
        ITEM_MORE = new GenomeListItem("More...", "", "More...");
    }

    /**
     * @param displayableName The name that can be shown to a user.
     * @param path            The location of the genome archive, can be a file path or URL
     * @param id              The id of the genome.
     */
    public GenomeListItem(String displayableName, String path, String id) {
        this.displayableName = displayableName;
        this.path = path;
        this.id = id;
    }

    public static GenomeListItem fromString(String str) {
        String[] tokens = str.split("\t");
        GenomeListItem item = new GenomeListItem(tokens[1], tokens[2], tokens[0]);
        item.hasDownloadedSequence = Boolean.parseBoolean(tokens[3]);
        return item;
    }


    public String printString() {
        return id + "\t" + displayableName + "\t" + path + "\t" + String.valueOf(hasDownloadedSequence);
    }

    public String getDisplayableName() {
        return displayableName;
    }

    public String getId() {
        return id;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public String getPath() {
        return path;
    }

    @Override
    public String toString() {
        return getDisplayableName();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GenomeListItem that = (GenomeListItem) o;

        if (displayableName != null ? !displayableName.equals(that.displayableName) : that.displayableName != null)
            return false;
        if (id != null ? !id.equals(that.id) : that.id != null) return false;
        if (path != null ? !path.equals(that.path) : that.path != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return Objects.hash(displayableName, path, id);
    }


}
