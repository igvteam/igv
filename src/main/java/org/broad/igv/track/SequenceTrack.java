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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.AminoAcidManager;
import org.broad.igv.feature.AminoAcidSequence;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Strand;
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

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static org.broad.igv.prefs.Constants.MAX_SEQUENCE_RESOLUTION;
import static org.broad.igv.prefs.Constants.SHOW_SEQUENCE_TRANSLATION;


/**
 * @author jrobinso
 */
public class SequenceTrack extends AbstractTrack implements IGVEventObserver {

    private static Logger log = Logger.getLogger(SequenceTrack.class);

    private static final int SEQUENCE_HEIGHT = 14;

    private static String NAME = "Sequence";

    private Map<String, LoadedDataInterval<SeqCache>> loadedIntervalCache = new HashMap(200);

    private Map<String, Boolean> sequenceVisible;

    private SequenceRenderer sequenceRenderer = new SequenceRenderer();

    //should translated aminoacids be shown below the sequence?
    private boolean shouldShowTranslation = true;

    Strand strand = Strand.POSITIVE;

    private Rectangle arrowRect;

    public SequenceTrack(String name) {
        super(null, name, name);
        setSortable(false);
        shouldShowTranslation = PreferencesManager.getPreferences().getAsBoolean(SHOW_SEQUENCE_TRANSLATION);
        loadedIntervalCache = Collections.synchronizedMap(new HashMap<>());
        sequenceVisible = Collections.synchronizedMap(new HashMap<>());
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

    public void receiveEvent(Object event) {

        if (event instanceof FrameManager.ChangeEvent) {

            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).getFrames();
            Map<String, LoadedDataInterval<SeqCache>> newCache = Collections.synchronizedMap(new HashMap<>());
            for (ReferenceFrame f : frames) {
                newCache.put(f.getName(), loadedIntervalCache.get(f.getName()));
            }
            loadedIntervalCache = newCache;


        } else {
            log.info("Unknown event type: " + event.getClass());
        }
    }

    private void refreshAminoAcids() {
         for(LoadedDataInterval<SeqCache> i : loadedIntervalCache.values()) {
             i.getFeatures().refreshAminoAcids();
         }
    }

    @Override
    public void renderName(Graphics2D graphics, Rectangle trackRectangle, Rectangle visibleRectangle) {

        Font font = FontManager.getFont(fontSize);

        boolean visible = this.sequenceVisible.values().stream().anyMatch(v -> v == true);

        if (visible) {
            graphics.setFont(font);
            int textBaseline = trackRectangle.y + 12;
            graphics.drawString(NAME, trackRectangle.x + 5, textBaseline);

            int rx = trackRectangle.x + trackRectangle.width - 20;
            arrowRect = new Rectangle(rx, trackRectangle.y + 2, 15, 10);
            drawArrow(graphics);

            //Show icon when translation non-standard
            if (AminoAcidManager.getInstance().getCodonTable().getId() != AminoAcidManager.STANDARD_TABLE_ID) {
                Font labFont = font.deriveFont(Font.BOLD);
                graphics.setFont(labFont);
                graphics.drawString("A", rx - 20, textBaseline);
                graphics.setFont(font);
            }

            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
        }
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
            return interval != null && interval.contains(frame);
        }
    }

    @Override
    public void load(ReferenceFrame referenceFrame) {

        String chr = referenceFrame.getChrName();
        final Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();

        Chromosome chromosome = currentGenome.getChromosome(chr);
        int start = (int) referenceFrame.getOrigin();
        final int chromosomeLength = chromosome.getLength();

        int end = (int) referenceFrame.getEnd();
        int w = end - start;

        // Expand a bit for panning and AA caluclation
        start = Math.max(0, start - w / 2 + 2);
        end = Math.min(end + w / 2 + 2, chromosomeLength);

        Genome genome = currentGenome;
        String sequence = new String(genome.getSequence(chr, start, end));

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
        cache.refreshAminoAcids();
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

        if (!sequenceVisible.containsKey(frameName)) sequenceVisible.put(frameName, false);  // Default value

        if (visible != sequenceVisible.get(frameName)) {
            sequenceVisible.put(frameName, visible);
            UIUtilities.invokeAndWaitOnEventThread(() -> {
                if (context.getPanel() != null) {
                    context.getPanel().revalidate();
                }
            });
        }
        if (visible) {
            LoadedDataInterval<SeqCache> sequenceInterval = loadedIntervalCache.get(frameName);
            if (sequenceInterval != null) {
                sequenceRenderer.setStrand(strand);
                sequenceRenderer.draw(sequenceInterval, context, rect, shouldShowTranslation, resolutionThreshold);
            }
        }
    }


    @Override
    public int getHeight() {
        boolean visible = this.sequenceVisible.values().stream().anyMatch(v -> v == true);
        return visible ? SEQUENCE_HEIGHT +
                (shouldShowTranslation ? SequenceRenderer.TranslatedSequenceDrawer.TOTAL_HEIGHT : 0) :
                0;
    }


    @Override
    public boolean handleDataClick(TrackClickEvent e) {
        setShouldShowTranslation(!shouldShowTranslation);
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

    public void setShouldShowTranslation(boolean shouldShowTranslation) {
        this.shouldShowTranslation = shouldShowTranslation;
        // Remember this choice
        PreferencesManager.getPreferences().put(SHOW_SEQUENCE_TRANSLATION, shouldShowTranslation);
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
        m2.setSelected(shouldShowTranslation);
        m2.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                setShouldShowTranslation(m2.isSelected());
                repaint();
                IGV.getInstance().clearSelections();
            }
        });

        menu.add(m1);
        menu.add(m2);

        final JMenu transTableMenu = new JMenu("Translation Table");
        for (AminoAcidManager.CodonTable codonTable : AminoAcidManager.getInstance().getAllCodonTables()) {
            JMenuItem item = getCodonTableMenuItem(codonTable);
            transTableMenu.add(item);
        }
        menu.add(transTableMenu);

        return menu;
    }

    private JCheckBoxMenuItem getCodonTableMenuItem(AminoAcidManager.CodonTable codonTable) {

        JCheckBoxMenuItem item = new JCheckBoxMenuItem();
        String fullName = codonTable.getDisplayName();
        String shortName = fullName;
        if (fullName.length() > 40) {
            shortName = fullName.substring(0, 37) + "...";
            item.setToolTipText(fullName);
        }
        item.setText(shortName);
        final AminoAcidManager.CodonTableKey curKey = codonTable.getKey();
        item.setSelected(curKey.equals(AminoAcidManager.getInstance().getCodonTable().getKey()));
        item.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                AminoAcidManager.getInstance().setCodonTable(curKey);
                SequenceTrack.this.refreshAminoAcids();
                repaint();
            }
        });
        return item;
    }


    private void repaint() {
        // TODO -- what's really needed is a repaint of all panels the sequence track intersects
        IGV.getMainFrame().repaint();
    }

    // SequenceTrack does not expose its renderer

    public Renderer getRenderer() {
        return null;
    }

    @Override
    public String getNameValueString(int y) {
        String nvs = "<html>" + super.getNameValueString(y);
        nvs += "<br>Translation Table: ";
        nvs += AminoAcidManager.getInstance().getCodonTable().getDisplayName();
        return nvs;
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
            this.posAA = posAA;
            this.negAA = negAA;
        }

        public void refreshAminoAcids() {
            int mod = start % 3;
            int n1 = normalize3(3 - mod);
            int n2 = normalize3(n1 + 1);
            int n3 = normalize3(n2 + 1);

            String sequence = new String(seq);
            AminoAcidSequence[] posAA = {
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n1, sequence.substring(n1)),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n2, sequence.substring(n2)),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.POSITIVE, start + n3, sequence.substring(n3))
            };
            this.posAA = posAA;

            final int len = sequence.length();
            AminoAcidSequence[] negAA = {
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n1)),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n2)),
                    AminoAcidManager.getInstance().getAminoAcidSequence(Strand.NEGATIVE, start, sequence.substring(0, len - n3))
            };
            this.negAA = negAA;
        }
    }
}
