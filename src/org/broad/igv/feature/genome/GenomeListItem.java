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

import org.broad.igv.util.HttpUtils;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipException;

/**
 * A container for specific genome information which can be used to
 * manage loaded genomes.
 */
public class GenomeListItem {

    private String displayableName;
    private String location;
    private String id;
    private Boolean hasDownloadedSequence = null;

    public static final GenomeListItem ITEM_MORE;

    static{
        ITEM_MORE = new GenomeListItem("More...", "", "More...");
    }

    /**
     *
     * @param displayableName The name that can be shown to a user.
     * @param location        The location of the genome archive, can be a file path or URL
     * @param id              The id of the genome.
     */
    public GenomeListItem(String displayableName, String location, String id) {
        this.displayableName = displayableName;
        this.location = location;
        this.id = id;
    }

    public String getDisplayableName() {
        return displayableName;
    }


    public String getId() {
        return id;
    }

    public String getLocation() {
        if(location == null){
            GenomeListItem newItem = GenomeManager.searchGenomeList(this.id, GenomeManager.getInstance().getServerGenomeArchiveList());
            if(newItem != null){
                this.displayableName = newItem.displayableName;
                this.location = newItem.location;
            }
        }
        return location;
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
        if (location != null ? !location.equals(that.location) : that.location != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = displayableName != null ? displayableName.hashCode() : 0;
        result = 31 * result + (location != null ? location.hashCode() : 0);
        result = 31 * result + (id != null ? id.hashCode() : 0);
        return result;
    }

    /**
     * Check if the genome being referred to points to a local (on this machine)
     * sequence, which was downloaded from a server. So a user-created genome
     * which points to a local fasta file will return false, but one created
     * by {@link GenomeManager#downloadWholeGenome(String, java.io.File, java.awt.Frame)}
     * will return true
     * @return
     */
    public boolean hasDownloadedSequence(){
        if(hasDownloadedSequence == null){
            try {
                hasDownloadedSequence = checkHasDownloadedSequence();
            } catch (IOException e) {
                e.printStackTrace();
                hasDownloadedSequence = false;
            }
        }
        return hasDownloadedSequence;
    }

    private boolean checkHasDownloadedSequence() throws IOException{
        if(this.location == null) return false;
        if(HttpUtils.isRemoteURL(this.location)) return false;

        if(FastaUtils.isFastaPath(this.location)){
            return !HttpUtils.isRemoteURL(this.location);
        }

        try {
            GenomeDescriptor descriptor = GenomeManager.parseGenomeArchiveFile(new File(this.location));
            return descriptor.hasCustomSequenceLocation() && !HttpUtils.isRemoteURL(descriptor.getSequenceLocation());
        } catch (ZipException e) {
            return false;
        }


    }
}
