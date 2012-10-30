package org.broad.igv.hic.tools;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.sql.*;
import java.util.regex.Pattern;

/**
 * @author jrobinso
 *         Date: 10/29/12
 *         Time: 10:50 AM
 */
public class HiCDBUtils {

    private static String DB_DRIVER;
    private static String DB_URL;
    private static String DB_USER;
    private static String DB_PASSWORD;

    public static void main(String[] args) throws IOException, SQLException {

        String cmd = args[0];
        if (cmd.equals("frag")) {
            String f = args[1];
            insertFragments(f);

        } else if (cmd.equals("annot")) {
            String f = args[1];
            insertAnnotationList(f);
            updateFragmentAnnotations();

        } else if (cmd.equals("update")) {
            updateFragmentAnnotations();

        } else {

        }
    }


    public static void insertAnnotationList(String annotListFile) throws IOException, SQLException {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(annotListFile);
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if (!nextLine.startsWith("#")) {
                    System.out.print("Processing " + nextLine);
                    insertAnnotations(nextLine);
                    System.out.println("   DONE");
                }
            }
        } finally {
            if (reader != null) reader.close();
        }
    }

    private static void insertAnnotations(String line) throws IOException, SQLException {

        Connection dbConnection = null;

        String insertAnnotationSql = "INSERT INTO IGV.ANNOTATION " +
                "(CHR, BEG, END, NAME, TYPE, SUBTYPE, CELL_TYPE, FILE_NAME, SCORE, ANTIBODY, SOURCE, REPLICATE_NUMBER, LINE) " +
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";

        PreparedStatement ps = null;

        String[] tk = Globals.whitespacePattern.split(line);
        String fragmentFile = tk[0];
        String celltype = tk[1];
        String altname = tk[2];
        String type = tk[3];
        String subtype = tk[4];
        String antibody = tk[5];
        String source = tk[6];
        String replicate = tk.length > 7 ? tk[7] : "";

        BufferedReader annotationReader = null;

        try {
            annotationReader = ParsingUtils.openBufferedReader(fragmentFile);

            dbConnection = getDBConnection();
            ps = dbConnection.prepareStatement(insertAnnotationSql);
            dbConnection.setAutoCommit(false);


            String nextLine;
            int count = 0;

            while ((nextLine = annotationReader.readLine()) != null) {

                String[] tokens = Globals.whitespacePattern.split(nextLine);
                String chr = tokens[0];
                int beg = Integer.parseInt(tokens[1]);
                int end = Integer.parseInt(tokens[2]);
                String name = tokens[3];
                int score = Integer.parseInt(tokens[4]);

                ps.setString(1, chr);
                ps.setInt(2, beg);
                ps.setInt(3, end);
                ps.setString(4, name);
                ps.setString(5, type);
                ps.setString(6, subtype);
                ps.setString(7, celltype);
                ps.setString(8, altname);
                ps.setInt(9, score);
                ps.setString(10, antibody);
                ps.setString(11, source);
                ps.setString(12, replicate);
                ps.setString(13, nextLine);
                ps.addBatch();

                count++;
                if (count % 1000 == 0) {
                    ps.executeBatch();
                    count = 0;
                }

            }

            if (count > 0) ps.executeBatch();

            dbConnection.commit();

            annotationReader.close();
            annotationReader = null;


        } finally {
            ps.close();
            dbConnection.close();
            if (annotationReader != null) annotationReader.close();
        }
    }


    public static void insertFragments(String fragmentFile) throws IOException, SQLException {

        Connection dbConnection = null;

        String insertTableSQL = "INSERT INTO IGV.FRAGMENT "
                + "(TYPE, CHR, BEG, END, IDX) VALUES"
                + "(?,?,?,?, ?)";
        PreparedStatement ps = null;


        BufferedReader fragmentReader = null;
        Pattern pattern = Pattern.compile("\\s");
        try {
            fragmentReader = new BufferedReader(new FileReader(fragmentFile));

            dbConnection = getDBConnection();
            ps = dbConnection.prepareStatement(insertTableSQL);
            dbConnection.setAutoCommit(false);


            String nextLine;
            while ((nextLine = fragmentReader.readLine()) != null) {
                String[] tokens = pattern.split(nextLine);


                // A hack, could use IGV's genome alias definitions
                String chr = getChrAlias(tokens[0]);
                System.out.println("Processing " + chr);
                int beg = 0;
                int idx = 0;
                for (int i = 1; i < tokens.length; i++) {
                    int end = Integer.parseInt(tokens[i]);

                    ps.setString(1, "HindIII");
                    ps.setString(2, chr);
                    ps.setInt(3, beg);
                    ps.setInt(4, end);
                    ps.setInt(5, idx);
                    ps.addBatch();

                    beg = end;
                    idx++;
                }

                ps.executeBatch();
            }
            dbConnection.commit();


        } finally {
            dbConnection.close();
            fragmentReader.close();
        }
    }


    public static void updateFragmentAnnotations() throws SQLException {

        String selectLastSQL = "SELECT MAX(ANNOTATION_ID) FROM FRAGMENT_ANNOTATION";
        String selectAnnotationSql = "SELECT ID, CHR, BEG, END FROM ANNOTATION WHERE ID > ?";
        String selectFragSql = "SELECT ID FROM FRAGMENT WHERE CHR = ? and BEG <= ? and END >= ?";
        String updateSQL = "INSERT INTO FRAGMENT_ANNOTATION (FRAGMENT_ID, ANNOTATION_ID) VALUES (?, ?)";

        Connection dbConnection = null;

        PreparedStatement lastIdPrepStat = null;
        PreparedStatement annotPrepStat = null;
        PreparedStatement fragPrepStat = null;
        PreparedStatement updatePrepStat = null;

        ResultSet lastIdRS = null;
        ResultSet annotRS = null;
        ResultSet fragRS = null;

        try {

            dbConnection = getDBConnection();
            dbConnection.setAutoCommit(false);
            annotPrepStat = dbConnection.prepareStatement(selectAnnotationSql);
            fragPrepStat = dbConnection.prepareStatement(selectFragSql);
            updatePrepStat = dbConnection.prepareStatement(updateSQL);
            int updateCount = 0;

            lastIdPrepStat = dbConnection.prepareStatement(selectLastSQL);
            lastIdRS = lastIdPrepStat.executeQuery();
            int lastIDProcessed = 0;
            if (lastIdRS.next()) lastIDProcessed = lastIdRS.getInt(1);
            System.out.println("Last id processed = " + lastIDProcessed);

            annotPrepStat.setInt(1, lastIDProcessed);

            annotRS = annotPrepStat.executeQuery();
            while (annotRS.next()) {

                int annotID = annotRS.getInt(1);
                String chr = annotRS.getString(2);
                int beg = annotRS.getInt(3);
                int end = annotRS.getInt(4);

                fragPrepStat.setString(1, chr);
                fragPrepStat.setInt(2, end);
                fragPrepStat.setInt(3, beg);
                fragRS = fragPrepStat.executeQuery();
                while (fragRS.next()) {
                    int fragID = fragRS.getInt(1);
                    updatePrepStat.setInt(1, fragID);
                    updatePrepStat.setInt(2, annotID);
                    updatePrepStat.addBatch();
                    updateCount++;

                    if (updateCount % 10 == 0) {
                        System.out.println("Updating");
                        updatePrepStat.executeBatch();
                        updateCount = 0;
                        dbConnection.commit();
                    }
                }
                fragRS.close();
                fragRS = null;
            }

            if (updateCount > 0) updatePrepStat.executeBatch();

            dbConnection.commit();

        } finally {
            annotRS.close();
            lastIdRS.close();
            if (fragRS != null) fragRS.close();

            lastIdPrepStat.close();
            annotPrepStat.close();
            updatePrepStat.close();
            fragPrepStat.close();
            dbConnection.close();
        }
    }


    private static Connection getDBConnection() {

        if (DB_DRIVER == null) {
            DB_DRIVER = System.getProperty("DB_DRIVER");
            DB_URL = System.getProperty("DB_URL");
            DB_USER = System.getProperty("DB_USER");
            DB_PASSWORD = System.getProperty("DB_PASSWORD");
        }

        Connection dbConnection = null;

        try {

            Class.forName(DB_DRIVER);

        } catch (ClassNotFoundException e) {

            System.out.println(e.getMessage());

        }

        try {

            dbConnection = DriverManager.getConnection(DB_URL, DB_USER, DB_PASSWORD);
            return dbConnection;

        } catch (SQLException e) {

            System.out.println(e.getMessage());

        }

        return dbConnection;

    }


    private static String getChrAlias(String token) {
        if (token.equals("MT")) {
            return "chrM";
        } else if (!token.startsWith("chr")) {
            return "chr" + token;
        } else {
            return token;
        }
    }


}
