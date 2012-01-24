package org.broad.igv.dev.affective;

import org.broad.igv.util.ParsingUtils;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Parser for Affective partipant files.
 * <p/>
 * Examples:
 * <p/>
 * Participant3.csv
 * ,Participant,Notes,Day Number,Date,More Relaxed than Usual,More Tense than Usual,Complained of Feeling Badly,Medication Change,Weekday
 * 8/2/11,3,0,1,3/22/11,Yes,No,No,No,Tuesday
 * 8/2/11,3,0,2,3/23/11,Yes,No,No,No,Wednesday
 * 8/3/11,3,0,3,3/24/11,No,No,No,No,Thursday
 * <p/>
 * ParticipantAnnotation3.csv
 * ,Participant,Staff,Date,Day Number,Onset Time,Offset Time,Type of Annotation,Mastered Activity?,Preferred Activity?,Description,Comments,DBR?,Antecedents,Perceived Reason(s) for Behavior
 * 8/2/11,3,M3,3/22/11,1,09:17:31,09:30:58,Activity,Yes,Yes,Bouncing on Ball,None,None,None,None
 * 8/2/11,3,M3,3/22/11,1,08:58:00,09:10:39,Activity,Yes,No,Soaking foot,None,None,None,None
 * 8/2/11,3,M3,3/22/11,1,09:11:01,09:12:56,Activity,Yes,No,Relaxation,None,None,None,None
 *
 * @author Jim Robinson
 * @date 1/22/12
 */
public class AffectiveAnnotationParser {

    static Pattern commaPattern = Pattern.compile(",");
    static Pattern colonPattern = Pattern.compile(":");
    static SimpleDateFormat dateFormat = new SimpleDateFormat("MM/dd/yy");

    int dateColumn = 4;
    int startTimeColumn = 6;
    int endTimeColumn = 7;
    int typeColumn = 8;
    int descriptionColumn = 10;

    String[] headings;


    /**
     * Return a map of description => list of annotations
     *
     * @param file
     * @return
     * @throws IOException
     */
    public Map<String, List<Annotation>> parse(String file) throws IOException {

        BufferedReader br = null;
        TreeMap<String, List<Annotation>> annotationMap = new TreeMap<String, List<Annotation>>();
        try {
            br = ParsingUtils.openBufferedReader(file);
            String headerLine = br.readLine();
            headings = commaPattern.split(headerLine);
            for (int i = 1; i < headings.length; i++) {
                if (headings[i].equals("Date")) {
                    dateColumn = i;
                } else if (headings[i].equals("Onset Time")) {
                    startTimeColumn = i;
                } else if (headings[i].equals("Offset Time")) {
                    endTimeColumn = i;
                } else if (headings[i].equals("Type of Annotation")) {
                    typeColumn = i;
                } else if (headings[i].equals("Description")) {
                    descriptionColumn = i;
                }
            }

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = commaPattern.split(nextLine);
                if (tokens.length != headings.length) {
                    System.out.println("Skipping: " + nextLine);
                    continue;
                }

                String dateString;
                try {
                    Date date = dateFormat.parse(tokens[dateColumn]);
                    dateString = AffectiveChromosome.dateFormat.format(date);
                } catch (ParseException e) {
                    System.out.println("Could not parse date: " + tokens[dateColumn]);
                    continue; // Skip linke
                }

                int startTime = toSeconds(tokens[startTimeColumn]);
                int endTime = toSeconds(tokens[endTimeColumn]);
                String type = tokens[typeColumn];
                String description = tokens[descriptionColumn];
                Annotation annotation = new Annotation(dateString, startTime, endTime, type, description,
                        headings, tokens);

                List<Annotation> annotList = annotationMap.get(type);
                if (annotList == null) {
                    annotList = new ArrayList<Annotation>();
                    annotationMap.put(type, annotList);
                }
                annotList.add(annotation);

            }
            return annotationMap;
        } finally {
            if (br != null) br.close();
        }
    }

    private static int toSeconds(String timeString) {
        String[] tmp = colonPattern.split(timeString);
        int h = Integer.parseInt(tmp[0]) - AffectiveUtils.START_TIME_HR;
        int m = Integer.parseInt(tmp[1]);
        int s = Integer.parseInt(tmp[2]);
        return ((h * 60) + m) * 60 + s;
    }



}
