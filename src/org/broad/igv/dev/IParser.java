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

import org.broad.igv.exceptions.ParserException;

/**
 * This interface acts as a bridge between different containers of data,
 * separating retrieving that data from processing. The motivation was integration
 * with SQL databases.
 * <p/>
 * Consider a tab-delimited text file. The standard method of parsing
 * is to read each line, splitting each line into a {@code String[]},
 * and reading each string into a float, double, string, or whatever.
 * {@link StringArrayParser#getString(TContainer, TIndex)}
 * has TContainer = String[], and TIndex = Integer, and simply
 * returns the String at location given by the Integer.
 * <p/>
 * Consider the same data stored in a SQL database. Each line stores
 * essentially the information, although we have the advantage that
 * SQL columns are typed. {@link org.broad.igv.dev.db.SQLLineParserByLabel}
 * reads ResultSet, and uses column names (String) to identify the location.
 * One could also use Integers to index the columns, if desired. This way
 * we retain type information on columns; the data is never treated as a String
 * (unless it is).
 * <p/>
 * User: jacob
 * Date: 2012-Aug-30
 *
 * @see org.broad.igv.dev.db.SQLLineParserByLabel
 * @see org.broad.igv.dev.StringArrayParser
 */
public interface IParser<TContainer, TIndex> {

    byte getByte(TContainer obj, TIndex index) throws ParserException;

    short getShort(TContainer obj, TIndex index) throws ParserException;

    int getInt(TContainer obj, TIndex index) throws ParserException;

    float getFloat(TContainer obj, TIndex index) throws ParserException;

    double getDouble(TContainer obj, TIndex index) throws ParserException;

    String getString(TContainer obj, TIndex index) throws ParserException;

    /**
     * Return the number of data objects contained.
     * For an array, this would be the length
     *
     * @param obj
     * @return
     * @throws ParserException
     */
    int size(TContainer obj) throws ParserException;
}
