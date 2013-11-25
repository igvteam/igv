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

import java.util.Map;

/**
 * Interface for a session persistable object.
 *
 * @author User: jrobinso
 * Date: Feb 23, 2010
 */
public interface Persistable {

    /**
     * Return object state as a map of key-value string pairs
     * @return
     */

    public Map<String, String> getPersistentState();

    /**
     * Restore object state from a map of key-value string pairs
     * @param values
     * @param version
     */
    public void restorePersistentState(Map<String, String> values, int version);
}
