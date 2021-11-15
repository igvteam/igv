package org.broad.igv.jbrowse;

import com.google.gson.Gson;
import htsjdk.tribble.Feature;
import org.broad.igv.bedpe.BedPEFeature;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantRenderer;

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
        Coord[] chords = new Coord[features.size()];
        int index = 0;
        for (BedPEFeature f : features) {
            chords[index++] = new Coord(f);
        }
        sendChordsToJBrowse(chords, trackName, color, "0.5");
    }

    public static void sendAlignmentsToJBrowse(List<Alignment> alignments, String trackName, Color color) {

        Coord[] chords = new Coord[alignments.size()];
        int index = 0;
        for (Alignment a : alignments) {
            chords[index++] = new Coord(a);
        }
        sendChordsToJBrowse(chords, trackName, color, "0.02");
    }

    public static void sendVariantsToJBrowse(List<Feature> variants, String trackName, Color color) {

        Coord[] chords = new Coord[variants.size()];
        int index = 0;
        for (Feature f : variants) {
            if (f instanceof Variant) {
                Variant v = (Variant) f;
                Map<String, Object> attrs = v.getAttributes();
                if (attrs.containsKey("CHR2") && attrs.containsKey("END")) {
                    chords[index++] = new Coord(v);
                }
            }
        }
        sendChordsToJBrowse(chords, trackName, color, "0.5");
    }

    public static void sendChordsToJBrowse(Coord[] chords, String trackName, Color color, String alpha) {

        String colorString = "rgba(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + "," + alpha + ")";
        CircViewTrack t = new CircViewTrack(chords, trackName, colorString);
        CircViewMessage message = new CircViewMessage("addChords", t);

        Gson gson = new Gson();
        String json = gson.toJson(message);
        SocketSender.send(json);
    }
}


class CircViewMessage {
    String message;
    Object data;
    public CircViewMessage(String message, Object data) {
        this.message = message;
        this.data = data;
    }
}

class CircViewTrack {
    String name;
    String color;
    Coord [] chords;

    public CircViewTrack(Coord[] chords, String name, String color) {
        this.name = name;
        this.color = color;
        this.chords = chords;
    }
}


class Mate {
    String refName;
    int start;
    int end;

    public Mate(String refName, int start, int end) {
        this.refName = refName;
        this.start = start;
        this.end = end;
    }
}

class Coord {
    String uniqueId;
    String color;
    String refName;
    int start;
    int end;
    Mate mate;

    public Coord(BedPEFeature f) {
        this.uniqueId = f.chr1 + ":" + f.start1 + "-" + f.end1 + "_" + f.chr2 + ":" + f.start2 + "-" + f.end2;
        this.refName = f.chr1.startsWith("chr") ? f.chr1.substring(3) : f.chr1;
        this.start = f.start1;
        this.end = f.end1;
        this.mate = new Mate(f.chr2.startsWith("chr") ? f.chr2.substring(3) : f.chr2, f.start2, f.end2);
        this.color = "rgba(0, 0, 255, 0.1)";
    }

    public Coord(Alignment a) {
        ReadMate mate = a.getMate();
        this.uniqueId = a.getReadName();
        this.refName = a.getChr().startsWith("chr") ? a.getChr().substring(3) : a.getChr();
        this.start = a.getStart();
        this.end = a.getEnd();
        this.mate = new Mate(mate.getChr().startsWith("chr") ? mate.getChr().substring(3) : mate.getChr(),
                mate.getStart(), mate.getStart() + 1);
        this.color = "rgba(0, 0, 255, 0.02)";
    }

    public Coord(Variant v) {

        Map<String, Object> attrs = v.getAttributes();
        String chr2 = shortName(attrs.get("CHR2").toString());
        int end2 = Integer.parseInt(attrs.get("END").toString());
        int start2 = end2 - 1;
        String chr1 = shortName(v.getChr());
        int start1 = v.getStart();
        int end1 = v.getEnd();

        this.uniqueId = chr1 + "_" + start1 + ":" + end1 + "-" + chr2 + "_" + start2 + ":" + end2;
        this.refName = chr1;
        this.start = start1;
        this.end = end1;
        this.mate = new Mate(chr2, start2, end2);
        this.color = "rgb(0,0,255)";
    }

    static String shortName(String chr) {
        return chr.startsWith("chr") ? chr.substring(3) : chr;
    }

}

class SocketSender {

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