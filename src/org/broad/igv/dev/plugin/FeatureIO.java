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

package org.broad.igv.dev.plugin;

import java.util.List;
import java.util.Map;

/**
 * Interface common for writing Features out and in
 * User: jacob
 * Date: 2012-Oct-18
 */
interface FeatureIO {

    /**
     * It may be the case that the output is processed differently, depending on the input.
     * We allow for that by providing the commands and arguments here. Implementations are not
     * required to do anything with the input
     *
     * @param commands    Command portions of the input, e.g. {"find", "."}
     * @param argumentMap Arguments with their values. e.g. "-name", "myFile"
     */
    void setInputs(List<String> commands, Map<Argument, Object> argumentMap);
}
