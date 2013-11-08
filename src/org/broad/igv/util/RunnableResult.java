/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util;

/**
 * Enum for indicating how a long running computation completed.
 *
 * @author jacob
 * @date 2013-Oct-23
 */
public enum RunnableResult {
    SUCCESS,
    FAILURE,
    CANCELLED;

    public boolean isSuccess() {
        return this == SUCCESS;
    }
}
