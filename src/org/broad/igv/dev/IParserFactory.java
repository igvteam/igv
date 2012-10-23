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

import org.broad.igv.dev.db.SQLLineParserByIndex;

import java.sql.ResultSet;

/**
 * User: jacob
 * Date: 2012-Oct-23
 */
public class IParserFactory {

    /**
     * Return an appropriate IParser based on the class of rs
     *
     * @param rs
     * @param <TContainer>
     * @return
     */
    public static <TContainer> IParser<TContainer, Integer> getIndexParser(TContainer rs) {
        IParser parser;
        if (rs instanceof ResultSet) {
            parser = new SQLLineParserByIndex();
        } else if (rs instanceof String[]) {
            parser = new StringArrayParser();
        } else {
            throw new IllegalArgumentException("Line must be a ResultSet or String[], not " + rs.getClass());
        }
        return parser;
    }
}
