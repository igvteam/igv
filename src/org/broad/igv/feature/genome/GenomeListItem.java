/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
