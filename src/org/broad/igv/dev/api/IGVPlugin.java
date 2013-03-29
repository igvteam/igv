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

package org.broad.igv.dev.api;

/**
 * Interface to implement to interact with IGV GUI.
 * @author  jacob
 * @date 2012-Dec-21
 * @api
 */

public interface IGVPlugin {

    /**
     * Called when the plugin is loaded. Should only
     * be called once during each IGV session
     */
    void init();
}
