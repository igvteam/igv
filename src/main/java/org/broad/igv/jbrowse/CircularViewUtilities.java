package org.broad.igv.jbrowse;

import com.google.gson.Gson;
import htsjdk.tribble.Feature;
import org.broad.igv.bedpe.BedPEFeature;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ChromosomeColors;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.vcf.MateVariant;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.List;
import java.util.Map;

public class CircularViewUtilities {

    public static boolean ping() {
        try {
            String response = SocketSender.send("{\"message\": \"ping\"}", true);
            return "OK".equals(response);
        } catch (Exception e) {
            return false;
        }
    }

    public static void sendBedpeToJBrowse(List<BedPEFeature> features, String trackName, Color color) {
        Chord[] chords = new Chord[features.size()];
        int index = 0;
        for (BedPEFeature f : features) {
            chords[index++] = new Chord(f);
        }
        sendChordsToJBrowse(chords, trackName, color, "0.5");
    }

    public static void sendAlignmentsToJBrowse(List<Alignment> alignments, String trackName, Color color) {

        Chord[] chords = new Chord[alignments.size()];
        int index = 0;
        for (Alignment a : alignments) {
            chords[index++] = new Chord(a);
        }
        sendChordsToJBrowse(chords, trackName, color, "0.02");
    }

    public static void sendVariantsToJBrowse(List<Feature> variants, String trackName, Color color) {

        Chord[] chords = new Chord[variants.size()];
        int index = 0;
        for (Feature f : variants) {
            if (f instanceof Variant) {
                Variant v = f instanceof MateVariant ? ((MateVariant) f).mate : (Variant) f;
                Map<String, Object> attrs = v.getAttributes();
                if (attrs.containsKey("CHR2") && attrs.containsKey("END")) {
                    chords[index++] = new Chord(v);
                }
            }
        }
        sendChordsToJBrowse(chords, trackName, color, "0.5");
    }

    public static void sendChordsToJBrowse(Chord[] chords, String trackName, Color color, String alpha) {

        // We can't know if an assembly has been set, or if it has its the correct one.
        changeGenome(GenomeManager.getInstance().getCurrentGenome());

        String colorString = "rgba(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + "," + alpha + ")";
        CircViewTrack t = new CircViewTrack(chords, trackName, colorString);
        CircViewMessage message = new CircViewMessage("addChords", t);

        Gson gson = new Gson();
        String json = gson.toJson(message);
        SocketSender.send(json);
    }

    public static void changeGenome(Genome genome) {
        List<String> wgChrNames = genome.getLongChromosomeNames();
        CircViewRegion[] regions = new CircViewRegion[wgChrNames.size()];
        int idx = 0;
        for (String chr : wgChrNames) {
            Chromosome c = genome.getChromosome(chr);
            int length = c.getLength();
            Color color = ChromosomeColors.getColor(chr);
            String colorString = "rgb(" + ColorUtilities.colorToString(color) + ")";
            regions[idx++] = new CircViewRegion(chr, length, colorString);
        }
        CircViewAssembly assm = new CircViewAssembly(genome.getId(), genome.getDisplayName(), regions);
        CircViewMessage message = new CircViewMessage("setAssembly", assm);

        Gson gson = new Gson();
        SocketSender.send(gson.toJson(message));

    }

    static class SocketSender {

        static String send(String json) {
            return send(json, false);
        }

        static String send(String json, boolean suppressErrors) {
            Socket socket = null;
            PrintWriter out = null;
            BufferedReader in = null;
            try {
                socket = new Socket("127.0.0.1", 1234);
                out = new PrintWriter(socket.getOutputStream(), true);
                in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

                out.println(json);
                out.flush();
                String response = in.readLine();
                return response;
            } catch (UnknownHostException e) {
                String err = "Unknown host exception: " + e.getMessage();
                if (!suppressErrors) System.err.println(err);
                return err;

            } catch (IOException e) {
                String message = "IO Exception: " + e.getMessage();
                if (!suppressErrors) System.err.println(message);
                return message;
            } finally {
                try {
                    in.close();
                    out.close();
                    socket.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}



/*


const MINIMUM_SV_LENGTH = 1000000;

        const circViewIsInstalled = () => CircularView.isInstalled();

        const shortChrName = (chrName) => {
        return chrName.startsWith("chr") ? chrName.substring(3) : chrName;
        }

        const makePairedAlignmentChords = (alignments, color) => {
        color = color || 'rgba(0, 0, 255, 0.02)'
        const chords = [];
        for (let a of alignments) {
        const mate = a.mate;
        if (mate && mate.chr && mate.position) {
        chords.push({
        uniqueId: a.readName,
        refName: shortChrName(a.chr),
        start: a.start,
        end: a.end,
        mate: {
        refName: shortChrName(mate.chr),
        start: mate.position - 1,
        end: mate.position,
        },
        color: color
        });
        }
        }
        return chords;
        }

        const makeBedPEChords = (features, color) => {

        color = color || 'rgb(0,0,255)';

        return features.map(v => {

        // If v is a whole-genome feature, get the true underlying variant.
        const f = v._f || v;

        return {
        uniqueId: `${f.chr1}:${f.start1}-${f.end1}_${f.chr2}:${f.start2}-${f.end2}`,
        refName: shortChrName(f.chr1),
        start: f.start1,
        end: f.end1,
        mate: {
        refName: shortChrName(f.chr2),
        start: f.start2,
        end: f.end2,
        },
        color: color,
        igvtype: 'bedpe'
        }
        })
        }


        const makeVCFChords = (features, color) => {

        color = color || 'rgb(0,0,255)';

        const svFeatures = features.filter(v => {
        const f = v._f || v;
        const isLargeEnough = f.info.CHR2 && f.info.END &&
        (f.info.CHR2 !== f.chr || Math.abs(Number.parseInt(f.info.END) - f.pos) > MINIMUM_SV_LENGTH);
        return isLargeEnough;
        });
        return svFeatures.map(v => {

        // If v is a whole-genome feature, get the true underlying variant.
        const f = v._f || v;

        const pos2 = Number.parseInt(f.info.END);
        const start2 = pos2 - 100;
        const end2 = pos2 + 100;

        return {
        uniqueId: `${f.chr}:${f.start}-${f.end}_${f.info.CHR2}:${f.info.END}`,
        refName: shortChrName(f.chr),
        start: f.start,
        end: f.end,
        mate: {
        refName: shortChrName(f.info.CHR2),
        start: start2,
        end: end2
        },
        color: color,
        igvtype: 'vcf'
        }
        })
        }

        const makeCircViewChromosomes = (genome) => {
        const regions = [];
        const colors = [];
        for (let chrName of genome.wgChromosomeNames) {
        const chr = genome.getChromosome(chrName);
        colors.push(getChrColor(chr.name));
        regions.push(
        {
        name: chr.name,
        bpLength: chr.bpLength
        }
        )
        }
        return regions;
        }


        function createCircularView(el, browser) {

        const circularView = new CircularView(el, {

        assembly: {
        name: browser.genome.id,
        id: browser.genome.id,
        chromosomes: makeCircViewChromosomes(browser.genome)
        },

        onChordClick: (feature, chordTrack, pluginManager) => {

        const f1 = feature.data;
        const f2 = f1.mate;
        const flanking = 2000;

        const l1 = new Locus({chr: browser.genome.getChromosomeName(f1.refName), start: f1.start, end: f1.end});
        const l2 = new Locus({chr: browser.genome.getChromosomeName(f2.refName), start: f2.start, end: f2.end});

        let loci;
        if ("alignment" === f1.igvtype) {   // append
        loci = this.currentLoci().map(str => Locus.fromLocusString(str));
        for (let l of [l1, l2]) {
        if (!loci.some(locus => {
        return locus.contains(l)
        })) {
        // add flanking
        l.start = Math.max(0, l.start - flanking);
        l.end += flanking;
        loci.push(l)
        }
        }
        } else {
        l1.start = Math.max(0, l1.start - flanking);
        l1.end += flanking;
        l2.start = Math.max(0, l2.start - flanking);
        l2.end += flanking;
        loci = [l1, l2];
        }

        const searchString = loci.map(l => l.getLocusString()).join(" ");
        browser.search(searchString);
        }
        });
        browser.circularView = circularView;
        circularView.hide();
        return circularView;
        }

        export {circViewIsInstalled, makeBedPEChords, makePairedAlignmentChords, makeVCFChords, createCircularView}


        */