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
import org.w3c.dom.*;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * Object representation of a single {@code table} element of
 * a database profile. Contains static method for parsing dbXML files
 * <p/>
 * User: jacob
 * Date: 2012-Oct-31
 */
public class DBTable {

    private static Logger log = Logger.getLogger(DBTable.class);

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
                String format = attr.getNamedItem("format").getTextContent();

                String posEndColName = Utilities.getNullSafe(attr, "posEndColName");
                String startColString = Utilities.getNullSafe(attr, "startColIndex");
                String endColString = Utilities.getNullSafe(attr, "endColIndex");
                String binColName = Utilities.getNullSafe(attr, "binColName");
                String baseQuery = Utilities.getNullSafe(attr, "baseQuery");

                int startColIndex = Integer.parseInt(startColString);
                int endColIndex = Integer.parseInt(endColString);

                //If present, retrieve list of columns
                Element tableElement = (Element) tableNode;

                NodeList columnNodes = tableElement.getElementsByTagName("column");
                Map<Integer, String> columnLabelMap = null;
                if (columnNodes.getLength() > 0) {
                    columnLabelMap = new HashMap<Integer, String>(columnNodes.getLength());
                    for (int col = 0; col < columnNodes.getLength(); col++) {
                        Node child = columnNodes.item(col);
                        NamedNodeMap colAttr = child.getAttributes();

                        int fileIndex = Integer.parseInt(colAttr.getNamedItem("fileIndex").getTextContent());
                        String colLabel = Utilities.getNullSafe(colAttr, "colLabel");
                        columnLabelMap.put(fileIndex, colLabel);

                    }
                }

                NodeList headerLineNodes = tableElement.getElementsByTagName("headerLine");
                List<String> headerLines = null;
                if (headerLineNodes.getLength() > 0) {
                    headerLines = new ArrayList<String>(headerLineNodes.getLength());
                    for (int hl = 0; hl < headerLineNodes.getLength(); hl++) {
                        Node child = headerLineNodes.item(hl);
                        headerLines.add(child.getTextContent().trim());
                    }
                }

                DBTable table = new DBTable(dbLocator, tableName, format, binColName, chromoColName, posStartColName,
                        posEndColName, startColIndex, endColIndex, columnLabelMap, baseQuery, headerLines);
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
            return createDBLocator(document);

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

    private final ResourceLocator dbLocator;
    private final String tableName;
    private final String format;
    private final String binColName;

    private final String chromoColName;
    private final String posStartColName;
    private final String posEndColName;
    private final int startColIndex;
    private final int endColIndex;

    private final Map<Integer, String> columnLabelMap;
    private final String baseQuery;
    private final List<String> headerLines;

    /**
     * Generally just intended for testing, where all we try
     * to do is get all data from a db and don't need anything fancy
     *
     * @param dbLocator
     * @param tableName
     */
    public static DBTable build(ResourceLocator dbLocator, String tableName) {
        return new DBTable(dbLocator, tableName, null, null, null, null, null, 1, Integer.MAX_VALUE - 1, null, null, null);
    }

    public DBTable(ResourceLocator dbLocator, String tableName, String format, String binColName,
                   String chromoColName, String posStartColName, String posEndColName, int startColIndex, int endColIndex,
                   Map<Integer, String> columnLabelMap, String baseQuery, List<String> headerLines) {
        this.dbLocator = dbLocator;
        this.tableName = tableName;
        this.format = format;
        this.binColName = binColName;
        this.chromoColName = chromoColName;
        this.posStartColName = posStartColName;
        this.posEndColName = posEndColName;
        this.startColIndex = startColIndex;
        this.endColIndex = endColIndex;
        this.columnLabelMap = columnLabelMap;
        this.baseQuery = baseQuery;
        this.headerLines = headerLines;
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

    public Map<Integer, String> getColumnLabelMap() {
        return columnLabelMap;
    }

    public List<String> getHeaderLines() {
        return headerLines;
    }

    /**
     * Return an array of column labels in specified ordinal positions
     *
     * @return
     */
    public static String[] columnMapToArray(Map<Integer, String> columnLabelMap) {
        List<Integer> arrayIndexes = new ArrayList<Integer>(columnLabelMap.keySet());
        Collections.sort(arrayIndexes);
        int minArrayIndex = arrayIndexes.get(0);
        int maxArrayIndex = arrayIndexes.get(arrayIndexes.size() - 1);
        int colCount = maxArrayIndex + 1;
        String[] tokens = new String[colCount];

        for (int cc = minArrayIndex; cc < maxArrayIndex; cc++) {
            tokens[cc] = columnLabelMap.get(cc);
        }
        return tokens;
    }
}
