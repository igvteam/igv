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

import org.apache.log4j.Logger;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXParseException;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Oct-31
 */
public class DBProfileReader {

    private static Logger log = Logger.getLogger(DBProfileReader.class);

    public static List<DBTable> parseProfile(String profilePath) {
        InputStream profileStream = null;
        String tableName = null;

        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            ResourceLocator dbLocator = createDBLocator(document);

            NodeList tableNodes = document.getElementsByTagName("table");

            List<DBTable> tableList = new ArrayList<DBTable>(tableNodes.getLength());

            for (int tnum = 0; tnum < tableNodes.getLength(); tnum++) {
                Node tableNode = tableNodes.item(tnum);
                NamedNodeMap attr = tableNode.getAttributes();

                tableName = attr.getNamedItem("name").getTextContent();
                String chromoColName = attr.getNamedItem("chromoColName").getTextContent();
                String posStartColName = attr.getNamedItem("posStartColName").getTextContent();
                String posEndColName = attr.getNamedItem("posEndColName").getTextContent();
                String format = attr.getNamedItem("format").getTextContent();
                String startColString = Utilities.getNullSafe(attr, "startColIndex");
                String endColString = Utilities.getNullSafe(attr, "endColIndex");
                String binColName = Utilities.getNullSafe(attr, "binColName");
                String baseQuery = Utilities.getNullSafe(attr, "baseQuery");

                int startColIndex = Integer.parseInt(startColString);
                int endColIndex = Integer.parseInt(endColString);

                //If present, retrieve list of columns
                NodeList columns = tableNode.getChildNodes();
                DBReader.ColumnMap columnMap = null;

                if (columns.getLength() > 0) {
                    columnMap = new DBReader.ColumnMap();

                    for (int col = 0; col < columns.getLength(); col++) {
                        Node column = columns.item(col);
                        NamedNodeMap colAttr = column.getAttributes();
                        //Whitespace gets in as child nodes
                        if (colAttr == null) continue;
                        int fileIndex = Integer.parseInt(colAttr.getNamedItem("fileIndex").getTextContent());

                        String colLabel = Utilities.getNullSafe(colAttr, "colLabel");
                        String colIndexStr = Utilities.getNullSafe(colAttr, "colIndex");

                        if (colLabel == null && colIndexStr == null) {
                            String msg = String.format("Error parsing column %d of %s", col, tableName);
                            msg += "colName and colIndex both null, at least 1 must be specified.";
                            throw new SAXParseException(msg, null);
                        } else if (colIndexStr != null) {
                            int colIndex = Integer.parseInt(colIndexStr);
                            columnMap.put(fileIndex, colIndex);
                        }
                        columnMap.put(fileIndex, colLabel);
                    }
                }

                DBTable table = new DBTable(dbLocator, tableName, format, binColName, chromoColName, posStartColName,
                        posEndColName, startColIndex, endColIndex, columnMap, baseQuery);
                tableList.add(table);
            }

            return tableList;

        } catch (Exception e) {
            String msg = "Error reading profile " + profilePath + ", table " + tableName;
            MessageUtils.showErrorMessage(msg, e);
            return null;
        } finally {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static ResourceLocator createDBLocator(Document document) {
        Node db = document.getElementsByTagName("database").item(0);
        NamedNodeMap attr = db.getAttributes();
        String host = attr.getNamedItem("host").getTextContent();
        String path = attr.getNamedItem("path").getTextContent();
        String subprotocol = attr.getNamedItem("subprotocol").getTextContent();

        String port = Utilities.getNullSafe(attr, "port");
        String username = Utilities.getNullSafe(attr, "username");
        String password = Utilities.getNullSafe(attr, "password");

        ResourceLocator locator = new ResourceLocator(DBManager.createConnectionURL(subprotocol, host, path, port));
        locator.setUsername(username);
        locator.setPassword(password);

        return locator;
    }

    /**
     * Creates a ResourceLocator from the dbxml file specified at
     * {@code profilePath}. Note that this locator is against the database,
     * and does not contain information on any particular table.
     *
     * @param profilePath
     * @return
     */
    static ResourceLocator createDBLocator(String profilePath) {
        InputStream profileStream = null;
        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            ResourceLocator locator = createDBLocator(document);
            return locator;

        } catch (Exception e) {
            log.error("Error creating DB Locator", e);
            return null;
        } finally {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException e) {
                log.error("Error closing profile stream", e);
            }
        }
    }

    /**
     * //     * Retrieve a reader from the XML profile located at {@code profilePath}.
     * //     * If {@code tableName == null}, all tables in the profile are loaded.
     * //     * TODO If {@code tableName == null}, the user is prompted to choose a table from the list
     * //     *
     * //     * @param profilePath
     * //     * @param tableName
     * //     * @return
     * //
     */
//    public static List<SQLCodecSource> getFromProfile(String profilePath, String tableName) {
//        List<DBTable> tableList = parseProfile(profilePath);
//        List<SQLCodecSource> sources = new ArrayList<SQLCodecSource>(tableList.size());
//
//        for (DBTable table : tableList) {
//            if (tableName == null || table.getTableName().equals(tableName)) {
//
//            }
//        }
//        return sources;
//    }


    public static class DBTable {
        private final ResourceLocator dbLocator;
        private final String tableName;
        private final String format;
        private final String binColName;

        private final String chromoColName;
        private final String posStartColName;
        private final String posEndColName;
        private final int startColIndex;
        private final int endColIndex;

        private final DBReader.ColumnMap columnMap;
        private final String baseQuery;

        public DBTable(ResourceLocator dbLocator, String tableName, String format, String binColName,
                       String chromoColName, String posStartColName, String posEndColName, int startColIndex, int endColIndex,
                       DBReader.ColumnMap columnMap, String baseQuery) {
            this.dbLocator = dbLocator;
            this.tableName = tableName;
            this.format = format;
            this.binColName = binColName;
            this.chromoColName = chromoColName;
            this.posStartColName = posStartColName;
            this.posEndColName = posEndColName;
            this.startColIndex = startColIndex;
            this.endColIndex = endColIndex;
            this.columnMap = columnMap;
            this.baseQuery = baseQuery;
        }

        public ResourceLocator getDbLocator() {
            return dbLocator;
        }

        public String getBinColName() {
            return binColName;
        }

        public String getChromoColName() {
            return chromoColName;
        }

        public DBReader.ColumnMap getColumnMap() {
            return columnMap;
        }

        public int getEndColIndex() {
            return endColIndex;
        }

        public String getFormat() {
            return format;
        }

        public String getPosEndColName() {
            return posEndColName;
        }

        public String getPosStartColName() {
            return posStartColName;
        }

        public int getStartColIndex() {
            return startColIndex;
        }

        public String getTableName() {
            return tableName;
        }

        public String getBaseQuery() {
            return baseQuery;
        }
    }
}
