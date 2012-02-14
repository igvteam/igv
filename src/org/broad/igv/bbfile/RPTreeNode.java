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

package org.broad.igv.bbfile;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Jan 14, 2010
 * Time: 11:16:03 AM
 * To change this template use File | Settings | File Templates.
 */
public interface RPTreeNode {

    // Identifies the node as a leaf node or a child (non-leaf) node.
    public boolean isLeaf();

    // Returns the chromosome bounds belonging to the entire node.
    public RPChromosomeRegion getChromosomeBounds();

     // Note: compareRegions returns the following values:
     //   -2 indicates chromosome region is completely below node region
     //   -1 indicates that chromosome region intersect node region from below
     //  0 means that chromosome region is inclusive to node region
     //  1 indicates chromosome region intersects node region from above
     //  2 indicates that this region is completely above that region
    public int compareRegions(RPChromosomeRegion chromosomeRegion);

    // Returns the number of items assigned to the node.
    public int getItemCount();

    // Returns the indexed node item.
    public  RPTreeNodeItem getItem(int index);

    // Inserts new node item according to bounds rank
    public boolean insertItem(RPTreeNodeItem item);

    // Deletes indexed node item
    public boolean deleteItem(int index);

    // prints the node items
    public void printItems();
}
