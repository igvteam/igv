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

package org.broad.igv.session;

/**
 * Interface for a session persistable object.  Implementations must be able to return object state as a map of
 * key-value string pairs, and restore state from such a map.
 *
 * // TODO -- this whole scheme could probably be more elegantly handled with annotations.
 *
 * @author User: jrobinso
 * Date: Feb 23, 2010
 */
public interface Persistable {

    /**
     * Return object state as a map of key-value string pairs
     * @return
     */

    public RecursiveAttributes getPersistentState();

    /**
     * Restore object state from a map of key-value string pairs
     * @param values
     */
    public void restorePersistentState(RecursiveAttributes values);
}
