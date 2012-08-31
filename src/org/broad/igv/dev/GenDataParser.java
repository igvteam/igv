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
