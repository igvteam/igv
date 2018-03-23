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

package org.broad.igv.dev;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * User: jacob
 * Date: 2012-Aug-30
 */
public class GenDataParser {

    public final byte getByte(ResultSet rs, int index) throws SQLException {
        return rs.getByte(index + 1);
    }

    public final byte getByte(ResultSet rs, String label) throws SQLException {
        return rs.getByte(label);
    }

    public final byte getByte(String[] array, int index) throws NumberFormatException {
        return Byte.valueOf(array[index]);
    }

    public final short getShort(ResultSet rs, int index) throws SQLException {
        return rs.getShort(index + 1);
    }

    public final short getShort(ResultSet rs, String label) throws SQLException {
        return rs.getShort(label);
    }

    public final short getShort(String[] array, int index) throws NumberFormatException {
        return Short.valueOf(array[index]);
    }

    public final int getInt(ResultSet rs, int index) throws SQLException {
        return rs.getInt(index + 1);
    }

    public final int getInt(ResultSet rs, String label) throws SQLException {
        return rs.getInt(label);
    }

    public final int getInt(String[] array, int index) throws NumberFormatException {
        return Integer.valueOf(array[index]);
    }

    public final double getDouble(ResultSet rs, int index) throws SQLException {
        return rs.getDouble(index + 1);
    }

    public final double getDouble(ResultSet rs, String label) throws SQLException {
        return rs.getDouble(label);
    }

    public final double getDouble(String[] array, int index) throws NumberFormatException {
        return Double.valueOf(array[index]);
    }

    public final float getFloat(ResultSet rs, int index) throws SQLException {
        return rs.getFloat(index + 1);
    }

    public final float getFloat(ResultSet rs, String label) throws SQLException {
        return rs.getFloat(label);
    }

    public final float getFloat(String[] array, int index) throws NumberFormatException {
        return Float.valueOf(array[index]);
    }

    public final String getString(ResultSet rs, int index) throws SQLException {
        return rs.getString(index + 1);
    }

    public final String getString(ResultSet rs, String label) throws SQLException {
        return rs.getString(label);
    }

    public final String getString(String[] array, int index) throws NumberFormatException {
        return array[index];
    }

//    public final byte getByte(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getByte(index + 1);
//        }else if(obj instanceof String[]){
//            return Byte.valueOf(((String[]) obj)[index]);
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
//    public final short getShort(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getShort(index + 1);
//        }else if(obj instanceof String[]){
//            return Short.valueOf(((String[]) obj)[index]);
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
//    public final int getInt(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getInt(index + 1);
//        }else if(obj instanceof String[]){
//            return Integer.valueOf(((String[]) obj)[index]);
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
//    public final float getFloat(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getFloat(index + 1);
//        }else if(obj instanceof String[]){
//            return Float.valueOf(((String[]) obj)[index]);
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
//    public final double getDouble(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getDouble(index + 1);
//        }else if(obj instanceof String[]){
//            return Double.valueOf(((String[]) obj)[index]);
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
//    public final String getString(Object obj, int index) throws SQLException{
//        if(obj instanceof ResultSet){
//            return ((ResultSet) obj).getString(index + 1);
//        }else if(obj instanceof String[]){
//            return ((String[]) obj)[index];
//        }else{
//            throw new IllegalArgumentException("Input must be a ResultSet or String[]");
//        }
//    }
//
}
