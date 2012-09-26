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

/**
 * A container for specific genome information which can be used to
 * manage loaded genomes.
 */
public class GenomeListItem {

    private String displayableName;
    private String location;
    private String id;

    /**
     * Constructor.
     *
     * @param displayableName The name that can be shown to a user.
     * @param url             The url of the genome archive.
     * @param id              The id of the genome.
     */
    public GenomeListItem(String displayableName, String url, String id) {

        this.displayableName = displayableName;
        this.location = url;
        this.id = id;
    }

    public String getDisplayableName() {
        return displayableName;
    }


    public String getId() {
        return id;
    }


    public String getLocation() {
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
}
