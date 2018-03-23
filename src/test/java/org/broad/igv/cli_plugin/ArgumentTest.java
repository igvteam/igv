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

package org.broad.igv.cli_plugin;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
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
        boolean eq = Utilities.objectEqual(a0.getType(), a1.getType()) &&
                Utilities.objectEqual(a0.getCmdArg(), a1.getCmdArg()) &&
                Utilities.objectEqual(a0.getDefaultValue(), a1.getDefaultValue()) &&
                Utilities.objectEqual(a0.getEncodingCodec(), a1.getEncodingCodec()) &&
                Utilities.objectEqual(a0.getId(), a1.getId()) &&
                Utilities.objectEqual(a0.getName(), a1.getName()) &&
                Utilities.objectEqual(a0.isOutput(), a1.isOutput());

        eq &= Arrays.deepEquals(a0.getLibPaths(), a1.getLibPaths());
        return eq;
    }


}
