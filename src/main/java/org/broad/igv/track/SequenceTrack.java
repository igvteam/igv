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


package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.event.IGVEvent;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.aa.AminoAcidManager;
import org.broad.igv.feature.aa.AminoAcidSequence;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.aa.CodonTable;
import org.broad.igv.feature.aa.CodonTableManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;


/**
 * @author jrobinso
 */
public class SequenceTrack extends AbstractTrack implements IGVEventObserver {

    private static Logger log = LogManager.getLogger(SequenceTrack.class);

    private static final int SEQUENCE_HEIGHT = 14;

    private static String NAME = "Sequence";

    private Map<String, LoadedDataInterval<SeqCache>> loadedIntervalCache = new HashMap(200);

    private SequenceRenderer sequenceRenderer = new SequenceRenderer();

    //should translated aminoacids be shown below the sequence?
    private boolean showTranslation = true;

    Strand strand = Strand.POSITIVE;

    private Rectangle arrowRect;

    public SequenceTrack(String name) {
        super(null, name, name);
        setSortable(false);
        showTranslation = PreferencesManager.getPreferences().getAsBoolean(SHOW_SEQUENCE_TRANSLATION);
        loadedIntervalCache = Collections.synchronizedMap(new HashMap<>());
        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
    }

    public static String getReverseComplement(String sequence) {
        char[] complement = new char[sequence.length()];
        int jj = complement.length;
        for (int ii = 0; ii < sequence.length(); ii++) {
            char c = sequence.charAt(ii);
            jj--;
            switch (c) {
                case 'T':
                    complement[jj] = 'A';
                    break;
                case 'A':
                    complement[jj] = 'T';
                    break;
                case 'C':
                    complement[jj] = 'G';
                    break;
                case 'G':
                    complement[jj] = 'C';
                    break;
                case 't':
                    complement[jj] = 'a';
                    break;
                case 'a':
                    complement[jj] = 't';
                    break;
                case 'c':
                    complement[jj] = 'g';
                    break;
                case 'g':
                    complement[jj] = 'c';
                    break;
                default:
                    complement[jj] = c;
            }
        }
        return new String(complement);
    }

    public void receiveEvent(IGVEvent event) {

        if (event instanceof FrameManager.ChangeEvent) {
            // Remove cache for discarded frames.  This seems a rather round-about way to do it.
            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).frames();
            Map<String, LoadedDataInterval<SeqCache>> newCache = Collections.synchronizedMap(new HashMap<>());
            for (ReferenceFrame f : frames) {
                newCache.put(f.getName(), loadedIntervalCache.get(f.getName()));
            }
            loadedIntervalCache = newCache;

        } else {
            log.warn("Unknown event type: " + event.getClass());
        }
    }

    private void refreshAminoAcids(CodonTable codonTable) {
        for (LoadedDataInterval<SeqCache> i : loadedIntervalCache.values()) {
            CodonTable intervalCodonTable = codonTable != null ?
                    codonTable :
                    CodonTableManager.getInstance().getCodonTableForChromosome(i.range.chr);
            if (intervalCodonTable != null) {
                i.getFeatures().refreshAminoAcids(intervalCodonTable);
            }
        }
    }

    @Override
    public void renderName(Graphics2D g, Rectangle trackRectangle, Rectangle visibleRectangle) {

        // Use local graphics -- this method corrupts graphics context when exporting to "png" files
        Graphics2D graphics = (Graphics2D) g.create();

        Font font = FontManager.getFont(getFontSize());

        boolean visible = isVisible();

        if (visible) {
            graphics.setFont(font);
            int textBaseline = trackRectangle.y + 12;
            graphics.drawString(NAME, trackRectangle.x + 5, textBaseline);

            int rx = trackRectangle.x + trackRectangle.width - 20;
            arrowRect = new Rectangle(rx, trackRectangle.y + 2, 15, 10);
            drawArrow(graphics);

            //Show icon when translation non-standard
//            if (AminoAcidManager.getInstance().getCodonTable().getId() != AminoAcidManager.STANDARD_TABLE_ID) {
//                Font labFont = font.deriveFont(Font.BOLD);
//                graphics.setFont(labFont);
//                graphics.drawString("A", rx - 20, textBaseline);
//                graphics.setFont(font);
//            }

            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
        }

        graphics.dispose();
    }

    private void drawArrow(Graphics2D graphics) {
        GraphicUtils.drawHorizontalArrow(graphics, arrowRect, strand == Strand.POSITIVE);
    }


    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {

        int resolutionThreshold = PreferencesManager.getPreferences().getAsInt(MAX_SEQUENCE_RESOLUTION);
        boolean visible = frame.getScale() < resolutionThreshold && !frame.getChrName().equals(Globals.CHR_ALL);

        if (!visible) {
            return true; // Nothing to paint
        } else {
            LoadedDataInterval<SeqCache> interval = loadedIntervalCache.get(frame.getName());
            boolean ready = interval != null && interval.contains(frame);
            return ready;
        }
    }

    @Override
    public void load(ReferenceFrame referenceFrame) {

        final Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
        String chr = currentGenome.getCanonicalChrName(referenceFrame.getChrName());

        int start = (int) referenceFrame.getOrigin();

        Chromosome chromosome = currentGenome.getChromosome(chr);
        if (chromosome == null) {
            return;
        }

        final int chromosomeLength = chromosome.getLength();

        int end = (int) referenceFrame.getEnd();
        int w = end - start;

        // Expand a bit for panning and AA caluclation
        start = Math.max(0, start - w / 2 + 2);
        end = Math.min(end + w / 2 + 2, chromosomeLength);

        Genome genome = currentGenome;
        byte [] seqBytes = genome.getSequence(chr, start, end);
        if(seqBytes == null) {
            return;

        }
        String sequence = new String(seqBytes);

        int mod = start % 3;
        int n1 = normalize3(3 - mod);
        int n2 = normalize3(n1 + 1);
        int n3 = normalize3(n2 + 1);

        // Now trim sequence to prevent dangling AAs
        int deltaStart = start == 0 ? 0 : 2;
        int deltaEnd = end == chromosomeLength ? 0 : 02;
        start += deltaStart;
        end -= deltaEnd;
        final int len = sequence.length();
        byte[] seq = sequence.substring(deltaStart, len - deltaEnd).getBytes();

        SeqCache cache = new SeqCache(start, seq);

        CodonTable codonTable = CodonTableManager.getInstance().getCodonTableForChromosome(chr);
        cache.refreshAminoAcids(codonTable);
        loadedIntervalCache.put(referenceFrame.getName(), new LoadedDataInterval<>(chr, start, end, cache));
    }

    private static int normalize3(int n) {
        return n == 3 ? 0 : n;
    }

    /**
     * Render the sequence, and optionally the 3 frame translation table
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect) {

        int resolutionThreshold = PreferencesManager.getPreferences().getAsInt(MAX_SEQUENCE_RESOLUTION);
        boolean visible = context.getReferenceFrame().getScale() < resolutionThreshold &&
                !context.getChr().equals(Globals.CHR_ALL);
        final String frameName = context.getReferenceFrame().getName();

        if (visible) {
            LoadedDataInterval<SeqCache> sequenceInterval = loadedIntervalCache.get(frameName);
            if (sequenceInterval != null) {
                sequenceRenderer.setStrand(strand);
                sequenceRenderer.draw(sequenceInterval, context, rect, showTranslation, resolutionThreshold);
            }
        }
    }

    @Override
    public boolean isVisible() {
        int resolutionThreshold = PreferencesManager.getPreferences().getAsInt(MAX_SEQUENCE_RESOLUTION);
        return FrameManager.getFrames().stream().anyMatch(frame -> (frame.getScale() < resolutionThreshold &&
                !frame.getChrName().equals(Globals.CHR_ALL)));
    }

    @Override
    public boolean isFilterable() {
        return false;
    }

    @Override
    public int getHeight() {
        return isVisible() ? SEQUENCE_HEIGHT +
                (showTranslation ? SequenceRenderer.TranslatedSequenceDrawer.TOTAL_HEIGHT : 0) :
                0;
    }


    @Override
    public boolean handleDataClick(TrackClickEvent e) {
        setShowTranslation(!showTranslation);
        Object source = e.getMouseEvent().getSource();
        if (source instanceof JComponent) {
            UIUtilities.invokeOnEventThread(() -> repaint());
        }
        return true;
    }

    @Override
    public void handleNameClick(final MouseEvent e) {
        if (arrowRect != null && arrowRect.contains(e.getPoint())) {
            flipStrand();
        }

    }

    private void flipStrand() {
        strand = (strand == Strand.POSITIVE ? Strand.NEGATIVE : Strand.POSITIVE);
        IGV.getInstance().clearSelections();
        repaint();

    }

    public void setStrand(Strand strandValue) {
        strand = strandValue;
        PreferencesManager.getPreferences().put(SEQUENCE_TRANSLATION_STRAND, strand.toString());
        IGV.getInstance().clearSelections();
        repaint();
    }

    public void setShowTranslation(boolean showTranslation) {
        this.showTranslation = showTranslation;
        PreferencesManager.getPreferences().put(SHOW_SEQUENCE_TRANSLATION, showTranslation);
        repaint();
    }


    /**
     * Override to return a specialized popup menu
     *
     * @return
     */
    @Override
    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        JMenuItem m1 = new JMenuItem("Flip strand");
        m1.addActionListener(e -> flipStrand());

        final JCheckBoxMenuItem m2 = new JCheckBoxMenuItem("Show translation");
        m2.setSelected(showTranslation);
        m2.addActionListener(e -> {
            setShowTranslation(m2.isSelected());
            IGV.getInstance().clearSelections();
        });

        menu.add(m1);
        menu.add(m2);

        final JMenu transTableMenu = new JMenu("Translation Table");
        transTableMenu.add(getCodonTableMenuItem(null));   // The "null" or default item
        for (CodonTable codonTable : CodonTableManager.getInstance().getAllCodonTables()) {
            JMenuItem item = getCodonTableMenuItem(codonTable);
            transTableMenu.add(item);
        }
        menu.add(transTableMenu);

        return menu;
    }

    private JCheckBoxMenuItem getCodonTableMenuItem(CodonTable codonTable) {

        JCheckBoxMenuItem item = new JCheckBoxMenuItem();

        if (codonTable == null) {
            item.setText("Default");
        } else {
            String fullName = codonTable.getDisplayName();
            String shortName = fullName;
            if (fullName.length() > 40) {
                shortName = fullName.substring(0, 37) + "...";
                item.setToolTipText(fullName);
            }
            item.setText(shortName);
        }

        final Integer selectedID = codonTable == null ? null : codonTable.getId();

        // Get the explicitly set codon table, if any
        CodonTable currentCodonTable = CodonTableManager.getInstance().getCurrentCodonTable();
        boolean selected = (selectedID == null && currentCodonTable == null) ||
                (selectedID != null && currentCodonTable != null && selectedID.equals(currentCodonTable.getId()));
        item.setSelected(selected);

        item.addActionListener(e -> {
            CodonTableManager.getInstance().setCurrentCodonTable(codonTable);
            SequenceTrack.this.refreshAminoAcids(codonTable);
            repaint();
        });
        return item;
    }

    // SequenceTrack does not expose its renderer

    public Renderer getRenderer() {
        return null;
    }

    @Override
    public String getTooltipText(int y) {
        CodonTable explicitlySelectedTable = CodonTableManager.getInstance().getCurrentCodonTable();
        if (explicitlySelectedTable != null) {
            String nvs = "<html>" + super.getTooltipText(y);
            nvs += "<br>Translation Table: ";
            nvs += explicitlySelectedTable.getDisplayName();
            return nvs;
        } else {
            return "";
        }
    }


    public Strand getStrand() {
        return this.strand;
    }


    public static class SeqCache {

        public int start;
        public byte[] seq;
        public AminoAcidSequence[] posAA;
        public AminoAcidSequence[] negAA;


        public SeqCache(int start, byte[] seq) {
            this.start = start;
            this.seq = seq;
        }

        /**
         * Refresh the amino acids with a specific codon table -- called from explicit menu choice.
         *
         * @param codonTable
         */
        public void refreshAminoAcids(CodonTable codonTable) {
            int n1 = start % 3;
            int n2 = (start + 1) % 3;
            int n3 = (start + 2) % 3;

            String sequence = new String(seq);
            AminoAcidSequence[] posAA = {
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n1, sequence.substring(n1), codonTable),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n2, sequence.substring(n2), codonTable),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n3, sequence.substring(n3), codonTable)
            };
            this.posAA = posAA;

            final int len = sequence.length();
            AminoAcidSequence[] negAA = {
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n1), codonTable),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n2), codonTable),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n3), codonTable)
            };
            this.negAA = negAA;
        }
    }

    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("shouldShowTranslation", String.valueOf(showTranslation));
        element.setAttribute("sequenceTranslationStrandValue", String.valueOf(strand));

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("shouldShowTranslation")) {
            this.showTranslation = Boolean.valueOf(element.getAttribute("shouldShowTranslation"));
        }
        if (element.hasAttribute("sequenceTranslationStrandValue")) {
            this.strand = Strand.fromString(element.getAttribute("sequenceTranslationStrandValue"));
        }
    }

}
