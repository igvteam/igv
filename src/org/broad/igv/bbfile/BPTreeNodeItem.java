
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

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Dec 23, 2009
 * Time: 11:43:43 AM
 * To change this template use File | Settings | File Templates.
 */
package org.broad.igv.bbfile;

import java.lang.Comparable;

/*
*   BPTreeNodeItem interface for storage of B+ tree node item information.
*
*   Note: The alpha-numeric key string is used for positional insertion of
*    node items and searching of the B+ tree.
* */

interface BPTreeNodeItem  {

    // Returns the child node item or leaf item index in the B+ tree.
    long getItemIndex();

    // Identifies the item as a leaf item or a child node item.
    boolean isLeafItem();

    // Returns key used to position the item in parent node item list.
    String getChromKey();

    // Returns true if keys match, returns false if keys do not match.
    boolean chromKeysMatch(String chromKey);

    // Prints the tree node items.
    void print();
}
