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

package org.broad.igv.util;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.Globals;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.util.HashSet;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Dec-26
 */
public class RuntimeUtilsTest extends AbstractHeadlessTest {

    @Test
    public void testSimpleRuntimeIO() throws Exception{
        String[] fullCmd = null;
        String expLine = null;
        int expLines = -1;
        if(Globals.IS_MAC || Globals.IS_LINUX){
            fullCmd = new String[]{"/bin/bash", "-c", "whereis bash"};
            expLine = "/bin/bash";
            expLines = 1;
        }else if(Globals.IS_WINDOWS){
            fullCmd = new String[]{"cmd","/C","ver"};
            expLine = "Microsoft Windows";
            expLines = 2;
        }

        Assume.assumeNotNull(fullCmd);
        Process process = RuntimeUtils.startExternalProcess(fullCmd, null, null);
        InputStream is = process.getInputStream();
        BufferedReader bis = new BufferedReader(new InputStreamReader(is));

        String line;
        boolean foundLine = false;
        int count = 0;

        process.waitFor();
        while((line = bis.readLine()) != null){
            System.out.println(line);
            count++;
            foundLine |= line.contains(expLine);
        }
        assertTrue(foundLine);
        assertEquals(expLines, count);
    }

    //Test our ability to write to stdin and read from stdout
    //XXX This test won't work on windows
    @Test
    public void testWriteStdIn() throws Exception {

        Assume.assumeTrue(!Globals.IS_WINDOWS);

        Runtime run = Runtime.getRuntime();
        Process pr = null;

        String msg = "Never have I ever \n thought twas ever thus \n a good movie \n thus I'm glad it's over";
        try {
            pr = run.exec("grep -i thus", null, null);
        } catch (IOException e) {
            e.printStackTrace();
        }

        BufferedReader in = new BufferedReader(new InputStreamReader(pr.getInputStream()));
        BufferedReader err = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(pr.getOutputStream())), true);
        out.println(msg);

        String line;

        out.flush();
        out.close();
        pr.waitFor();

        int readlines = 0;
        while ((line = in.readLine()) != null) {
            //System.out.println(line);
            readlines++;
        }

        in.close();


        //System.out.println("errors:");

        int errlines = 0;
        while ((line = err.readLine()) != null) {
            System.out.println(line);
            errlines++;
        }

        assertEquals(2, readlines);
        assertEquals(0, errlines);
    }

    //Behavior is platform dependent
    @Ignore
    @Test
    public void testGetSize_String() throws Exception{
        String testObj = "abcdefghi";
        //Strings have some overhead: pointer itself, int offset,hash,count
        //Not exactly sure how much memory a char[] takes up, doesn't seem precisely linear
        long expSize = 32 + 32*3 + 32 + testObj.length()*16 + 8;

        long actSize = JavaAgent.getObjectSizeRecursive(testObj, new HashSet<Object>());
        assertEquals(String.format("Characters: %d. Difference: %d in size", testObj.length(), expSize - actSize), expSize, actSize);
    }
}
