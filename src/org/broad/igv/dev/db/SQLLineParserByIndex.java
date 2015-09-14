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
