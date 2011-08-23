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

package org.broad.igv.feature.genome;

/**
 * A container for specific genome information which can be used to
 * manage loaded genomes.
 */
public class GenomeListItem {

    private String displayableName;
    private String location;
    private String id;
    private boolean userDefined = false;

    /**
     * Constructor.
     *
     * @param displayableName The name that can be shown to a user.
     * @param url             The url of the genome archive.
     * @param id              The id of the genome.
     * @param isUserDefined
     */
    public GenomeListItem(String displayableName, String url, String id, boolean isUserDefined) {

        this.displayableName = displayableName;
        this.location = url;
        this.id = id;
        this.userDefined = isUserDefined;
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


    public boolean isUserDefined() {
        return userDefined;
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

        if (userDefined != that.userDefined) return false;
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
        result = 31 * result + (userDefined ? 1 : 0);
        return result;
    }
}
