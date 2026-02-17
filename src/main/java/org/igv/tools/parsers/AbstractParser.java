
package org.igv.tools.parsers;

import org.igv.track.DataType;

import java.util.Hashtable;

/**
 * @author jrobinso
 */
public abstract class AbstractParser {

    private static Hashtable<String, String> chromosomeNames = new Hashtable();

    private DataConsumer dataConsumer;
    String trackLine = null;
    private DataType dataType = DataType.COPY_NUMBER;
    private String[] headings;

    public AbstractParser(DataConsumer dataConsumer) {
        this.dataConsumer = dataConsumer;
    }

    /**
     * @return the dataConsumer
     */
    public DataConsumer getDataConsumer() {
        return dataConsumer;
    }

    protected void setTrackParameters() {
        dataConsumer.setTrackParameters(dataType, trackLine, headings);
    }

    /**
     * Note:  This is an exact copy of the method in ExpressionFileParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     */
    void parseComment(String comment) {
        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            trackLine = tmp;
        } else {
            String[] tokens = tmp.split("=");
            if (tokens.length != 2) {
                return;
            }
            String key = tokens[0].trim().toLowerCase();
            if (key.equals("type")) {
                try {
                    setTrackType(DataType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {
                }
            }
        }
    }

    /**
     * @param headings the headings to set
     */
    public void setHeadings(String[] headings) {
        this.headings = headings;
    }

    /**
     * @return the headings
     */
    public String[] getHeadings() {
        return headings;
    }

    /**
     * @return the trackType
     */
    public DataType getTrackType() {
        return dataType;
    }

    /**
     * @param dataType the trackType to set
     */
    public void setTrackType(DataType dataType) {
        this.dataType = dataType;
    }

}
