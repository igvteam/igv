package org.igv.tools;

import java.io.*;
import java.util.List;

/**
 * Sum the values from a list of wig files and output a new wig file with the totals.  The wig files have to have
 * exactly the same coordinates row per row.
 * <p/>
 * TODO:  NOTE -- this only works with variable step at the moment.  The first line in each file must be variableStep
 *
 * @author jrobinso
 * @date Mar 16, 2011
 */
public class WigSummer {

    public static void sumWigs(List<File> inputs, File output) throws IOException {

        PrintWriter out = null;
        BufferedReader[] in = new BufferedReader[inputs.size()];

        try {
            out = new PrintWriter(new BufferedWriter(new FileWriter(output)));
            for (int i = 0; i < inputs.size(); i++) {
                in[i] = new BufferedReader(new FileReader(inputs.get(i)));
            }

            int lineNumber = 1;
            String firstLine = null;
            for (BufferedReader reader : in) {
                if (firstLine == null) {
                    firstLine = reader.readLine().trim();
                    if (!firstLine.startsWith("variableStep")) {
                        throw new RuntimeException("First line must be a variableStep line");
                    }
                } else {
                    String tmp = reader.readLine().trim();
                    if (!tmp.equals(firstLine)) {
                        throw new RuntimeException("First line of all input files must be equal");
                    }
                }

                out.println(firstLine);
                lineNumber++;
            }

            String tmp = null;
            while ((tmp = in[0].readLine()) != null) {

                String[] tokens = tmp.split("\t");
                int position = Integer.parseInt(tokens[0]);
                float v = Float.parseFloat(tokens[1]);

                for (int i = 1; i < in.length; i++) {
                    tmp = in[i].readLine();
                    tokens = tmp.split("\t");
                    int p = Integer.parseInt(tokens[0]);

                    if (p != position) {
                        throw new RuntimeException("Positions not all equal at line number " + lineNumber);
                    }

                    v += Float.parseFloat(tokens[1]);
                }

                out.println(position + "\t" + v);

            }
        } finally {
            out.close();
            for(BufferedReader reader : in) {
                reader.close();
            }
        }


    }
}
