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

package org.broad.igv.batch;

import org.apache.commons.lang.StringUtils;
import org.broad.igv.dev.plugin.batch.Command;

import java.util.List;

/**
 * User: jacob
 * Date: 2012-Nov-30
 */
public class EchoCommand implements Command {
    @Override
    public String run(List<String> args) {
        return StringUtils.join(args, " ");
    }
}
