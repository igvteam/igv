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

package org.broad.igv.dev.db;

import org.broad.igv.dev.IParser;
import org.broad.igv.exceptions.ParserException;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * User: jacob
 * Date: 2012-Aug-30
 *
 * @see IParser
 * @see org.broad.igv.dev.StringArrayParser
 */
public class SQLLineParserByLabel implements IParser<ResultSet, String> {

    @Override
    public final byte getByte(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getByte(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final short getShort(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getShort(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final int getInt(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getInt(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final double getDouble(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getDouble(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final float getFloat(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getFloat(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final String getString(ResultSet rs, String label) throws ParserException {
        try {
            return rs.getString(label);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public int size(ResultSet obj) throws ParserException {
        try {
            return obj.getMetaData().getColumnCount();
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

}
