
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
