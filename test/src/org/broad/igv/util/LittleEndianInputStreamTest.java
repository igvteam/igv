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

import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;

import org.junit.Ignore;


/**
 * @author jrobinso
 *         Date: 7/19/12
 *         Time: 8:54 AM
 */
@Ignore
public class LittleEndianInputStreamTest {


    @BeforeClass
    public static void setup() throws Exception {
        createTestFile();
    }

    @Test
    public void testRead() throws Exception {
        LittleEndianInputStream lis = new LittleEndianInputStream(new BufferedInputStream(new FileInputStream("les_test.bin")));
        assertEquals("Binary test file", lis.readString());
        assertEquals(Float.MAX_VALUE, lis.readFloat());
        assertEquals(Byte.MAX_VALUE, lis.readByte());
        assertEquals(Short.MAX_VALUE, lis.readShort());
        assertEquals(Integer.MAX_VALUE, lis.readInt());
        assertEquals(Long.MAX_VALUE, lis.readLong());
        assertEquals(Double.MAX_VALUE, lis.readDouble());
        lis.close();
    }


    static void createTestFile() throws IOException {

        LittleEndianOutputStream los = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream("les_test.bin")));

        los.writeString("Binary test file");
        los.writeFloat(Float.MAX_VALUE);
        los.writeByte(Byte.MAX_VALUE);
        los.writeShort(Short.MAX_VALUE);
        los.writeInt(Integer.MAX_VALUE);
        los.writeLong(Long.MAX_VALUE);
        los.writeDouble(Double.MAX_VALUE);

        los.close();


    }
}
