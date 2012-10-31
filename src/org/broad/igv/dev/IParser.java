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
