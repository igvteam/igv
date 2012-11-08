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
 * Retrieve value from a SQL line. Note that all input {@code index}
 * values should be 0-based, NOT the SQL standard 1-based
 * User: jacob
 * Date: 2012-Aug-30
 *
 * @see IParser
 * @see org.broad.igv.dev.StringArrayParser
 */
public class SQLLineParserByIndex implements IParser<ResultSet, Integer> {

    @Override
    public final byte getByte(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getByte(index + 1);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final short getShort(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getShort(index + 1);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final int getInt(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getInt(index + 1);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final double getDouble(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getDouble(index + 1);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final float getFloat(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getFloat(index + 1);
        } catch (SQLException e) {
            throw new ParserException(e.getMessage(), -1);
        }
    }

    @Override
    public final String getString(ResultSet rs, Integer index) throws ParserException {
        try {
            return rs.getString(index + 1);
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
