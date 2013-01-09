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

package org.broad.igv.cli_plugin;

import com.google.java.contract.util.Objects;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-27
 */

public class ArgumentTest extends AbstractHeadlessTest{

    @Test
    public void testTypeText() throws Exception{
        File loclib = new File("test/lib/RuntimeUtils.jar");
        String[] libURLs = new String[]{"file://" + loclib.getAbsolutePath(), loclib.getPath(), "http://www.example.com/test.jar"};
        Argument inArg = new Argument("name", Argument.InputType.TEXT, "cmdArg", "defVal", "encCodec", libURLs, true, "id");

        Argument outArg = TestUtils.marshallUnmarshall(inArg);
        assertTrue(argumentsEqual(inArg, outArg));
    }

    private boolean argumentsEqual(Argument a0, Argument a1) throws Exception{
        boolean eq = Objects.equal(a0.getType(), a1.getType()) &&
                Objects.equal(a0.getCmdArg(), a1.getCmdArg()) &&
                Objects.equal(a0.getDefaultValue(), a1.getDefaultValue()) &&
                Objects.equal(a0.getEncodingCodec(), a1.getEncodingCodec()) &&
                Objects.equal(a0.getId(), a1.getId()) &&
                Objects.equal(a0.getName(), a1.getName()) &&
                Objects.equal(a0.isOutput(), a1.isOutput());

        eq &= Arrays.deepEquals(a0.getLibPaths(), a1.getLibPaths());
        return eq;
    }


}
