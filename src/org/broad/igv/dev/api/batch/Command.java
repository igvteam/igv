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

package org.broad.igv.dev.api.batch;

import java.util.List;

/**
 * Interface which must be implemented by
 * any class which is a custom batch script command
 * <p/>
 * User: jacob
 * Date: 2012-Nov-30
 * @api
 */

public interface Command {

    /**
     * @param args The string arguments passed AFTER the command name. May be empty. <br/>
     *             Example:<br/>
     *             myCommand abba zabba
     *             <br/>
     *             would provide args = ['abba', 'zabba']
     *             <br/>
     *             myCommand
     *             <br/>
     *             would provide an empty list
     * @return A status string, which should be "OK" if everything went fine, and start with "ERROR" if something went wrong.
     *         Never null
     */
    String run(List<String> args);
}
