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

package org.broad.igv.dev;

import org.broad.igv.exceptions.ParserException;

/**
 * User: jacob
 * Date: 2012-Aug-30
 */
public class StringArrayParser implements IParser<String[], Integer> {

    @Override
    public final byte getByte(String[] array, Integer index) throws NumberFormatException {
        return Byte.parseByte(array[index].trim());
    }

    @Override
    public final short getShort(String[] array, Integer index) throws NumberFormatException {
        return Short.parseShort(array[index].trim());
    }

    @Override
    public final int getInt(String[] array, Integer index) throws NumberFormatException {
        return (int) Double.parseDouble(array[index].trim());
    }

    @Override
    public final double getDouble(String[] array, Integer index) throws NumberFormatException {
        return Double.parseDouble(array[index].trim());
    }

    @Override
    public final float getFloat(String[] array, Integer index) throws NumberFormatException {
        return Float.parseFloat(array[index].trim());
    }

    @Override
    public final String getString(String[] array, Integer index) {
        return array[index];
    }

    @Override
    public int size(String[] obj) throws ParserException {
        return obj.length;
    }

}
