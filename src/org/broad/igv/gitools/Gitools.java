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

package org.broad.igv.gitools;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManagerUI;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.RegionScoreType;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.collections.DoubleArrayList;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.Socket;
import java.util.*;

/**
 * Class for integrating with Gitools. Writes out data in TDM format.
 * We implement GeneListUIActionListener so the gene list can be selected using
 * our current UI
 * @author jrobinso
 *         Date: 11/13/12
 *         Time: 9:26 PM
 */
public class Gitools implements IGVPlugin{

    private static Logger log = Logger.getLogger(Gitools.class);

    public static final String ENABLE_PROPERTY = "enableGitools";

    private static int port = 50151;
    private static String host = "localhost";

    @Override
    public void init() {
        addMenuItems();
    }

    private static void addMenuItems(){
        boolean showTDMButton = Boolean.parseBoolean(System.getProperty(Gitools.ENABLE_PROPERTY, "true"));
        if(showTDMButton){
            JMenu gitoolsMenu = new JMenu("Gitools Heatmaps");

            JMenuItem directLoadItem = new JMenuItem("Load Gene Matrix in Gitools");
            directLoadItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {

                    if (Gitools.canConnect() ) {
                        GeneListManagerUI dialog = GeneListManagerUI.getInstance(IGV.getMainFrame(),
                                "Gitools Load", new Gitools.DirectLoadListener());
                        dialog.setVisible(true);
                    } else {
                        JOptionPane.showMessageDialog(IGV.getMainFrame(), "To be able to browse the gene matrix you need to install and open Gitools.\n Download it from http://www.gitools.org.");
                    }

                }
            });
            gitoolsMenu.add(directLoadItem);

            JMenuItem gitoolsItem = new JMenuItem("Export Gene Matrix (TDM)...");
            gitoolsItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    GeneListManagerUI dialog = GeneListManagerUI.getInstance(IGV.getMainFrame(),
                            "Export TDM", new Gitools.ExportFileListener());
                    dialog.setVisible(true);
                }
            });
            gitoolsMenu.add(gitoolsItem);
            IGV.getInstance().addOtherToolMenu(gitoolsMenu);
        }
    }

    public static void main(String[] args){
        try {
            sendCommand("version");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Send a command to gitools, and read the response.
     * @param command
     * @return First line of response only TODO This is stupid but readLine is hanging so hack it for now
     * @throws IOException
     */
    private static List<String> sendCommand(String command) throws IOException{
        Socket socket = null;
        try{
            socket = new Socket(host, port);
            socket.setSoTimeout(1000);
            PrintWriter out = new PrintWriter(socket.getOutputStream(),
                    true);
            BufferedReader in = new BufferedReader(new InputStreamReader(
                    socket.getInputStream()));

            out.println(command);
            List<String> response = new ArrayList<String>();
            String line;
            while( (line = in.readLine()) != null){
                response.add(line);
                break;
            }
            return response;

        }catch(IOException e){
            System.out.println(e);
            throw new IOException("Error communicating with gitools", e);
        }finally{
            if(socket != null) socket.close();
        }
    }

    public static boolean canConnect(){
        boolean canConnect = false;
        try{
            List<String> version = sendCommand("version");
            canConnect = version.get(0).toLowerCase().contains("gitools");
        } catch (IOException e) {
            //
        }
        return canConnect;
    }

    /**
     * Export data at specified loci to temp file, and tell gitools to load it.
     * gitools must already be running
     * @param lociStrings
     * @return
     * @throws IOException
     */
    public static List<String> gitoolsLoad(String name, List<String> lociStrings) throws IOException{
        String prefix = name + "-igv";
        prefix = prefix.replace(" ", "_");
        File tmpFile = File.createTempFile(prefix, ".tdm");

        exportTDM(lociStrings, tmpFile);

        return sendCommand(String.format("load %s", tmpFile.getAbsolutePath()));
    }

    public static void exportTDM(List<String> lociStrings, File file) throws IOException {

        // Convert the loci strings to a list of loci, if the loci represents multiple features (e.g. isoforms) use the largest
        int averageFeatureSize = 0;
        List<NamedFeature> loci = new ArrayList<NamedFeature>(lociStrings.size());
        for (String l : lociStrings) {
            NamedFeature feature = FeatureDB.getFeature(l);
            if (feature == null) {
                feature = Locus.fromString(l);
            }
            if (feature != null) {
                loci.add(feature);
                averageFeatureSize += (feature.getEnd() - feature.getStart());
            }
        }
        if (loci.size() > 0) averageFeatureSize /= loci.size();

        // Determine data types -- all data tracks + mutation, and samples
        LinkedHashSet<TrackType> loadedTypes = new LinkedHashSet<TrackType>();

        List<Track> tracks = IGV.getInstance().getAllTracks();
        for (Track t : tracks) {
            if ((t instanceof DataTrack || t.getTrackType() == TrackType.MUTATION) && t.isVisible()) {
                loadedTypes.add(t.getTrackType());
                //samples.add(t.getSample());
            }
        }

        // set an appropriate zoom level for this feature set.  Very approximate
        int zoom = 0;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null) {
            double averageChrLength = genome.getTotalLength() / genome.getChromosomes().size();
            zoom = (int) (Math.log(averageChrLength / averageFeatureSize) / Globals.log2) + 1;
        }


        // Loop though tracks and loci and gather data by sample & gene
        Map<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        for (Track t : tracks) {
            if (!t.isVisible()) continue;

            String sampleName = t.getSample();
            List<Track> overlays = IGV.getInstance().getOverlayTracks(t);

            for (NamedFeature locus : loci) {

                double regionScore;
                if (t instanceof DataTrack) {
                    DataTrack dataTrack = (DataTrack) t;
                    regionScore = dataTrack.getAverageScore(locus.getChr(), locus.getStart(), locus.getEnd(), zoom);
                    addToSampleData(sampleDataMap, sampleName, locus.getName(), t.getTrackType(), regionScore);
                }

                if(overlays != null) {
                    for (Track overlay : overlays) {
                        if (overlay.getTrackType() == TrackType.MUTATION) {
                            regionScore = overlay.getRegionScore(locus.getChr(), locus.getStart(), locus.getEnd(), zoom,
                                    RegionScoreType.MUTATION_COUNT, locus.getName());
                            //Only add if we found a mutation. Should we put it in anyway?
                            if (regionScore > 0) {
                                addToSampleData(sampleDataMap, sampleName, locus.getName(), overlay.getTrackType(), regionScore);
                            }
                        }
                    }
                }
            }

        }

        writeTDM(loadedTypes, sampleDataMap, file);
    }

    /**
     * Export calculated data to file
     * @param loadedTypes Ordered set of track types
     * @param sampleDataMap Map from sample name to SampleData
     * @param file Output file
     * @throws IOException
     */
    private static void writeTDM(LinkedHashSet<TrackType> loadedTypes, Map<String, SampleData> sampleDataMap, File file) throws IOException {
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));

            pw.print("Sample\tLocus");
            for (TrackType tt : loadedTypes) {
                pw.print("\t" + tt.name());
            }
            pw.println();

            for (SampleData sd : sampleDataMap.values()) {
                pw.print(sd.sample + "\t" + sd.locus);
                for (TrackType tt : loadedTypes) {
                    pw.print('\t');
                    double[] values = sd.getValues(tt);
                    if (values == null) {
                        pw.print("-");
                    } else {
                        double avg = StatUtils.max(values);
                        pw.print(avg);
                    }
                }
                pw.println();
            }
        } catch(IOException e){
            throw new IOException("Error exporting TDM", e);
        } finally{
            if (pw != null) pw.close();
        }
    }

    private static void addToSampleData(Map<String, SampleData> sampleDataMap, String sampleName, String locusString, TrackType tt, double regionScore){
        if (!Double.isNaN(regionScore)) {
            //String locusString = locus.getName();
            String key = sampleName + "_" + locusString;

            SampleData sd = sampleDataMap.get(key);
            if (sd == null) {

                sd = new SampleData(sampleName, locusString);
                sampleDataMap.put(key, sd);
            }
            sd.addValue(tt, regionScore);
        }
    }

    static class SampleData {

        String sample;
        String locus;
        Map<TrackType, DoubleArrayList> valueMap;

        SampleData(String sample, String locus) {
            this.sample = sample;
            this.locus = locus;
            this.valueMap = new HashMap<TrackType, DoubleArrayList>();
        }

        void addValue(TrackType tt, double value) {

            DoubleArrayList vlist = valueMap.get(tt);
            if (vlist == null) {
                vlist = new DoubleArrayList();
                valueMap.put(tt, vlist);
            }
            vlist.add(value);
        }

        public double[] getValues(TrackType tt) {
            DoubleArrayList dal = valueMap.get(tt);
            return dal == null ? null : dal.toArray();
        }
    }

    public static class ExportFileListener implements GeneListManagerUI.GeneListListener{

        @Override
        public void actionPerformed(JDialog dialog, GeneList geneList) {
            File file = FileDialogUtils.chooseFile("Export TDM file", null, FileDialogUtils.SAVE);

            // Check file TDM extension
            String currentExtension = FilenameUtils.getExtension(file.getName());
            if (!currentExtension.equalsIgnoreCase("TDM")) {
                file = new File(file.getAbsolutePath() + ".tdm");
            }

            if (file != null) {
                try {
                    Gitools.exportTDM(geneList.getLoci(), file);
                } catch (IOException exc) {
                    MessageUtils.showErrorMessage("Error exporting TDM", exc);
                }
            }
        }
    }

    public static class DirectLoadListener implements GeneListManagerUI.GeneListListener{

        @Override
        public void actionPerformed(JDialog dialog, GeneList geneList) {

            //Test communication
            boolean canConnect = canConnect();
            if(!canConnect){
                MessageUtils.showMessage("Cannot communicate with gitools, check that it is running on port " + port);
                return;
            }

            try {
                gitoolsLoad(geneList.getName(), geneList.getLoci());
            } catch (IOException exc) {
                MessageUtils.showErrorMessage(exc.getMessage(), exc);
            }


        }
    }
}
