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

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.sql.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Class for reading data from SQL database, where the lines
 * are of a format we can read with existing codecs
 *
 * @author Jacob Silterra
 * @date 29 May 2012
 */
public class SQLCodecReader extends DBReader {

    private static Logger log = Logger.getLogger(SQLCodecReader.class);

    protected String template;
    protected FeatureCodec codec;

    public SQLCodecReader(String template, FeatureCodec codec){
       this.codec = codec;


//        this.template = template;
//        InputStream templateStream = SQLCodecReader.class.getResourceAsStream("resources/" + template);
//        try {
//            Document document = Utilities.createDOMDocumentFromXmlStream(templateStream);
//            //Since we are parsing a stored resource, these should only happen if the file
//            //gets corrupted
//        } catch (ParserConfigurationException e) {
//            e.printStackTrace();
//        } catch (SAXException e) {
//            e.printStackTrace();
//        } catch (IOException e) {
//            e.printStackTrace();
//        } finally {
//            try {
//                templateStream.close();
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }

        }


    public Iterable<Feature> load(ResourceLocator locator, Genome genome) {

        List<Feature> featureList;

        ResultSet rs = null;
        Statement st = null;
        Connection conn = null;
        try {
            conn = DBManager.getConnection(locator.getUrl());
            String query = locator.getDescription();
            st = conn.createStatement();
            rs = st.executeQuery(query);

            featureList = new ArrayList<Feature>(rs.getFetchSize());
            while (rs.next()) {
                Feature feat = this.parseLine(rs);
                featureList.add(feat);
            }


            System.out.println("Disconnected from database");
        } catch (SQLException e) {
            log.error("Database error", e);
            throw new RuntimeException("Database error", e);
        } finally {
            closeResources(rs, st, conn);
        }

        return featureList;
    }

    /**
     * Turn a line from a result set into a feature
     * TODO Make abstract
     * @param rs
     * @return
     * @throws SQLException
     */
    protected Feature parseLine(ResultSet rs) throws SQLException{
        String[] tokens = lineToArray(rs);
        //TODO GET RID OF THIS, IT'S BAD AND I FEEL BAD FOR WRITING IT -JS
        String line = StringUtils.join(tokens, "\t");
        return codec.decode(line);

    }

    /**
     * Convert a the current line to an array of strings
     * @param rs
     * @return
     * @throws SQLException
     */
    protected String[] lineToArray(ResultSet rs) throws SQLException {
        int colCount = rs.getMetaData().getColumnCount();
        String[] tokens = new String[colCount];
        for(int cc= 0; cc < colCount; cc++){
            //SQL indexes from 1
            tokens[cc] = rs.getString(cc+1);
        }
        return tokens;
    }

}
