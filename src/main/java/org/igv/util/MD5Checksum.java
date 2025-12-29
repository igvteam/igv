package org.igv.util;

/**
 * Created by IntelliJ IDEA.
 * User: jesse
 * Date: Jan 19, 2010
 * Time: 5:30:36 PM
 * To change this template use File | Settings | File Templates.
 */

import org.igv.exceptions.DataLoadException;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.security.MessageDigest;

public class MD5Checksum {

    private static byte[] createChecksum(String filename) throws
            Exception {
        InputStream fis = new FileInputStream(filename);

        byte[] buffer = new byte[1024];
        MessageDigest complete = MessageDigest.getInstance("MD5");
        int numRead;
        do {
            numRead = fis.read(buffer);
            if (numRead > 0) {
                complete.update(buffer, 0, numRead);
            }
        } while (numRead != -1);
        fis.close();
        return complete.digest();
    }

    public static String getMD5Checksum(String filename) throws Exception {
        File fileName = new File(filename);

        //TODO Make a better output file naming convention.
        if (fileName.isFile()) {
            byte[] b = createChecksum(filename);
            String result = "";
            for (int i = 0; i < b.length; i++) {
                result +=
                        Integer.toString((b[i] & 0xff) + 0x100, 16).substring(1);
            }
            return result;
        }
        throw new DataLoadException("No file found", filename);
    }
}
