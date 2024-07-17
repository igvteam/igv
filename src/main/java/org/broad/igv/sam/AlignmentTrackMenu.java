package org.broad.igv.sam;

import htsjdk.samtools.SAMTag;
import org.broad.igv.Globals;
import org.broad.igv.event.AlignmentTrackEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.Strand;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.mods.BaseModficationFilter;
import org.broad.igv.sam.mods.BaseModificationKey;
import org.broad.igv.sam.mods.BaseModificationUtils;
import org.broad.igv.sashimi.SashimiPlot;
import org.broad.igv.tools.PFMExporter;
import org.broad.igv.track.SequenceTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.InsertSizeSettingsDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.supdiagram.SupplementaryAlignmentDiagramDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.extview.ExtendViewClient;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import static org.broad.igv.prefs.Constants.*;

/**
 * Popup menu class for AlignmentTrack.  The menu gets instantiated from TrackPanelComponent on right-click in the
 * alignment track or its associated name panel.
 */
class AlignmentTrackMenu extends IGVPopupMenu {

    private static final Logger log = LogManager.getLogger(AlignmentTrackMenu.class);

    private final AlignmentTrack alignmentTrack;
    private final AlignmentDataManager dataManager;
    private final AlignmentTrack.RenderOptions renderOptions;
    private static int nClusters = 2;

    AlignmentTrackMenu(AlignmentTrack alignmentTrack, final TrackClickEvent e) {

        this.alignmentTrack = alignmentTrack;
        this.dataManager = alignmentTrack.getDataManager();
        this.renderOptions = alignmentTrack.getRenderOptions();
        final Alignment clickedAlignment = alignmentTrack.getAlignmentAt(e);

        // Title
        JLabel popupTitle = new JLabel("  " + alignmentTrack.getName(), JLabel.CENTER);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);

        // Circular view items -- optional
        if (CircularViewUtilities.ping()) {
            addSeparator();
            JMenuItem item = new JMenuItem("Add Discordant Pairs to Circular View");
            item.setEnabled(alignmentTrack.getDataManager().isPairedEnd());
            add(item);
            item.addActionListener(ae -> alignmentTrack.sendPairsToCircularView(e));

            JMenuItem item2 = new JMenuItem("Add Split Reads to Circular View");
            add(item2);
            item2.addActionListener(ae -> alignmentTrack.sendSplitToCircularView(e));
        }

        // Some generic items from TrackMenuUtils
        Collection<Track> tracks = List.of(alignmentTrack);
        addSeparator();
        add(TrackMenuUtils.getTrackRenameItem(tracks));

        JMenuItem item = new JMenuItem("Change Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(tracks));
        add(item);

        // Experiment type  (RNA, THIRD GEN, OTHER)
        addSeparator();
        addExperimentTypeMenuItem();

        // Group, sort, color, shade, and pack
        addSeparator();
        addGroupMenuItem(e);
        addSortMenuItem();
        addColorByMenuItem();
        addShadeAlignmentsMenuItem();
        //addFilterMenuItem();
        addPackMenuItem();

        // Shading and mismatch items
        addSeparator();
        addShadeBaseByMenuItem();
        JMenuItem misMatchesItem = addShowMismatchesMenuItem();
        JMenuItem showAllItem = addShowAllBasesMenuItem();
        misMatchesItem.addActionListener(new Deselector(misMatchesItem, showAllItem));
        showAllItem.addActionListener(new Deselector(showAllItem, misMatchesItem));

        // Duplicates
        addDuplicatesMenuItem();

        // Hide small indels
        JMenuItem smallIndelsItem = new JCheckBoxMenuItem("Hide small indels");
        smallIndelsItem.setSelected(renderOptions.isHideSmallIndels());
        smallIndelsItem.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            if (smallIndelsItem.isSelected()) {
                String sith = MessageUtils.showInputDialog("Small indel threshold: ", String.valueOf(renderOptions.getSmallIndelThreshold()));
                if (sith == null) {
                    // dialogue was cancelled so no change should be made
                    return;
                } else {
                    try {
                        renderOptions.setSmallIndelThreshold(Integer.parseInt(sith));
                    } catch (NumberFormatException exc) {
                        log.error("Error setting small indel threshold - not an integer", exc);
                        //error so no change should be made
                        return;
                    }
                }
            }
            renderOptions.setHideSmallIndels(smallIndelsItem.isSelected());
            alignmentTrack.repaint();
        }));
        add(smallIndelsItem);

        // Paired end items
        if (dataManager.isPairedEnd()) {
            addSeparator();
            addViewAsPairsMenuItem();
            if (clickedAlignment != null) {
                addGoToMate(e, clickedAlignment);
                showMateRegion(e, clickedAlignment);
            }
            addInsertSizeMenuItem();
        }

        // Third gen (primarily) items
        addSeparator();
        addThirdGenItems(clickedAlignment, e);

        if(AlignmentTrack.ExperimentType.SBX == alignmentTrack.getExperimentType()) {
            addSBXItems(clickedAlignment, e);
        }

        // Display mode items
        addSeparator();
        TrackMenuUtils.addDisplayModeItems(tracks, this);

        // Select alignment items
        addSeparator();
        addSelectByNameItem(alignmentTrack, e);
        addClearSelectionsMenuItem();

        // Copy items
        addSeparator();
        addCopyToClipboardItem(e, clickedAlignment);
        addCopySequenceItems(e);
        addConsensusSequence(e);

        // Blat items
        addSeparator();
        addBlatItem(e);
        addBlatClippingItems(e);

        // Insertion items, only if clicked over an insertion
        AlignmentBlock insertion = alignmentTrack.getInsertion(clickedAlignment, e.getMouseEvent().getX());
        if (insertion != null) {
            addSeparator();
            addInsertionItems(insertion);
        }

        // Sashimi plot
        if (alignmentTrack.getExperimentType() == AlignmentTrack.ExperimentType.RNA) {
            addSeparator();
            JMenuItem sashimi = new JMenuItem("Sashimi Plot");
            sashimi.addActionListener(e1 -> SashimiPlot.openSashimiPlot());
            sashimi.setEnabled(alignmentTrack.getExperimentType() == AlignmentTrack.ExperimentType.RNA);
            add(sashimi);
        }

        // Experimental items
        if (alignmentTrack.getExperimentType() == AlignmentTrack.ExperimentType.THIRD_GEN) {
            addSeparator();
            addClusterItem(e);
        }

        // Show alignments, coverage, splice junctions
        addSeparator();
        addShowItems();
    }

    private void addDuplicatesMenuItem() {
        JMenu duplicatesMenu = new JMenu("Duplicates");
        for (AlignmentTrack.DuplicatesOption option : AlignmentTrack.DuplicatesOption.values()) {
            JRadioButtonMenuItem mi = new JRadioButtonMenuItem(option.label);
            final AlignmentTrack.DuplicatesOption previous = renderOptions.getDuplicatesOption();
            mi.setSelected(previous == option);
            mi.addActionListener(aEvt -> {
                renderOptions.setDuplicatesOption(option);
                if (previous != option) {
                    if (previous.filtered != option.filtered) {
                        // duplicates are filtered out when loading the read data so a reload has to be performed in this case
                        IGVEventBus.getInstance().post(new AlignmentTrackEvent(AlignmentTrackEvent.Type.RELOAD));
                    } else {
                        alignmentTrack.repaint();
                    }
                } else {
                    alignmentTrack.repaint();
                }
            });

            duplicatesMenu.add(mi);
        }
        add(duplicatesMenu);
    }

    private void addShowChimericRegions(final AlignmentTrack alignmentTrack, final TrackClickEvent e, final Alignment clickedAlignment) {

        JMenuItem item = new JMenuItem("View chimeric alignments in split screen");
        if (clickedAlignment != null && clickedAlignment.getAttribute(SAMTag.SA.name()) != null) {
            item.setEnabled(true);
            item.addActionListener(aEvt -> {
                final String saTag = clickedAlignment.getAttribute(SAMTag.SA.name()).toString();
                try {
                    List<SupplementaryAlignment> supplementaryAlignments = SupplementaryAlignment.parseFromSATag(saTag);
                    alignmentTrack.setSelectedAlignment(clickedAlignment);
                    FrameManager.addNewLociToFrames(e.getFrame(), supplementaryAlignments, alignmentTrack.getSelectedReadNames().keySet());
                } catch (final Exception ex) {
                    MessageUtils.showMessage("Failed to handle SA tag: " + saTag + " due to " + ex.getMessage());
                    item.setEnabled(false);
                }
            });
        } else {
            item.setEnabled(false);
        }
        add(item);
    }

    private void addShowDiagram(final TrackClickEvent e, final Alignment clickedAlignment) {
        JMenuItem item = new JMenuItem("Supplementary Reads Diagram");
        if (clickedAlignment != null && clickedAlignment.getAttribute(SAMTag.SA.name()) != null) {
            item.setEnabled(true);
            item.addActionListener(aEvt -> {
                try {
                    final SupplementaryAlignmentDiagramDialog frame = new SupplementaryAlignmentDiagramDialog(IGV.getInstance().getMainFrame(), clickedAlignment, new Dimension(500, 250));
                    frame.setVisible(true);
                } catch (final Exception ex) {
                    MessageUtils.showMessage("Failed to handle SA tag: " + clickedAlignment.getAttribute(SAMTag.SA.name()) + " due to " + ex.getMessage());
                    item.setEnabled(false);
                }
            });
        } else {
            item.setEnabled(false);
        }
        add(item);
    }


    private void addClusterItem(TrackClickEvent e) {

        JMenuItem item = new JMenuItem("Cluster alignments  *EXPERIMENTAL*");

        final ReferenceFrame frame;
        if (e.getFrame() == null && FrameManager.getFrames().size() == 1) {
            frame = FrameManager.getFrames().get(0);
        } else {
            frame = e.getFrame();
        }

        item.setEnabled(frame != null);
        add(item);

        item.addActionListener(ae -> {
            //This shouldn't ever be true, but just in case it's more user-friendly
            if (frame == null) {
                MessageUtils.showMessage("Unknown region bounds");
                return;
            }

            String nString = MessageUtils.showInputDialog("Enter the number of clusters", String.valueOf(nClusters));
            if (nString == null) {
                return;
            }
            try {
                nClusters = Integer.parseInt(nString);
            } catch (NumberFormatException e1) {
                MessageUtils.showMessage("Clusters size must be an integer");
                return;
            }

            final int start = (int) frame.getOrigin();
            final int end = (int) frame.getEnd();

            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            ClusterUtils clusterUtils = new ClusterUtils(interval);
            boolean success = clusterUtils.clusterAlignments(frame.getChrName(), start, end, nClusters);

            if (success) {
                groupAlignments(AlignmentTrack.GroupOption.CLUSTER, null, null);
                alignmentTrack.repaint();
            }

        });


    }


    /**
     * Item for exporting "consensus" sequence of region, based on loaded alignments.
     *
     * @param e
     */
    private void addConsensusSequence(TrackClickEvent e) {

        JMenuItem item = new JMenuItem("Copy consensus sequence");

        final ReferenceFrame frame;
        if (e.getFrame() == null && FrameManager.getFrames().size() == 1) {
            frame = FrameManager.getFrames().get(0);
        } else {
            frame = e.getFrame();
        }

        item.setEnabled(frame != null);
        add(item);

        item.addActionListener(ae -> {

            if (frame == null) {  // Should never happen
                MessageUtils.showMessage("Unknown region bounds, cannot export consensus");
                return;
            }
            final int start = (int) frame.getOrigin();
            final int end = (int) frame.getEnd();
            if ((end - start) > 1000000) {
                MessageUtils.showMessage("Cannot export region more than 1 Megabase");
                return;
            }
            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            AlignmentCounts counts = interval.getCounts();
            String text = PFMExporter.createPFMText(counts, frame.getChrName(), start, end);
            StringUtils.copyTextToClipboard(text);
        });


    }

    private JMenu getBisulfiteContextMenuItem(ButtonGroup group) {

        JMenu bisulfiteContextMenu = new JMenu("bisulfite mode");

        JRadioButtonMenuItem nomeESeqOption = null;
        boolean showNomeESeq = alignmentTrack.getPreferences().getAsBoolean(SAM_NOMESEQ_ENABLED);
        if (showNomeESeq) {
            nomeESeqOption = new JRadioButtonMenuItem("NOMe-seq bisulfite mode");
            nomeESeqOption.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.NOMESEQ);
            nomeESeqOption.addActionListener(aEvt -> {
                alignmentTrack.setColorOption(AlignmentTrack.ColorOption.NOMESEQ);
                alignmentTrack.repaint();
            });
            group.add(nomeESeqOption);
        }

        for (final AlignmentTrack.BisulfiteContext item : AlignmentTrack.BisulfiteContext.values()) {

            JRadioButtonMenuItem m1 = new JRadioButtonMenuItem(item.getLabel());
            m1.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.BISULFITE && renderOptions.bisulfiteContext == item);
            m1.addActionListener(aEvt -> {
                alignmentTrack.setColorOption(AlignmentTrack.ColorOption.BISULFITE);
                alignmentTrack.setBisulfiteContext(item);
                alignmentTrack.repaint();
            });
            bisulfiteContextMenu.add(m1);
            group.add(m1);
        }

        if (nomeESeqOption != null) {
            bisulfiteContextMenu.add(nomeESeqOption);
        }

        return bisulfiteContextMenu;

    }

    void addSelectByNameItem(AlignmentTrack track, TrackClickEvent e) {
        JMenuItem item = new JMenuItem("Select by name...");
        final Alignment alignment = track.getAlignmentAt(e);
        String alignmentName = alignment == null ? "" : alignment.getReadName();
        item.addActionListener(aEvt -> {
            String val = MessageUtils.showInputDialog("Enter read name: ", alignmentName);
            if (val != null && val.trim().length() > 0) {
                alignmentTrack.getSelectedReadNames().put(val, alignmentTrack.getReadNamePalette().get(val));
                alignmentTrack.repaint();
            }
        });
        add(item);
    }

    void addExperimentTypeMenuItem() {
        Map<String, AlignmentTrack.ExperimentType> mappings = new LinkedHashMap<>();
        mappings.put("Other", AlignmentTrack.ExperimentType.OTHER);
        mappings.put("RNA", AlignmentTrack.ExperimentType.RNA);
        mappings.put("3rd Gen", AlignmentTrack.ExperimentType.THIRD_GEN);
        mappings.put("SBX", AlignmentTrack.ExperimentType.SBX);
        //mappings.put("Bisulfite", ExperimentType.BISULFITE);
        JMenu groupMenu = new JMenu("Experiment Type");
        ButtonGroup group = new ButtonGroup();
        for (Map.Entry<String, AlignmentTrack.ExperimentType> el : mappings.entrySet()) {
            JCheckBoxMenuItem mi = getExperimentTypeMenuItem(el.getKey(), el.getValue());
            groupMenu.add(mi);
            group.add(mi);
        }
        add(groupMenu);
    }

    private JCheckBoxMenuItem getExperimentTypeMenuItem(String label, final AlignmentTrack.ExperimentType option) {
        JCheckBoxMenuItem mi = new JCheckBoxMenuItem(label);
        mi.setSelected(alignmentTrack.getExperimentType() == option);
        mi.addActionListener(aEvt -> alignmentTrack.setExperimentType(option));
        return mi;
    }

    void addGroupMenuItem(final TrackClickEvent te) {//ReferenceFrame frame) {
        final MouseEvent me = te.getMouseEvent();
        ReferenceFrame frame = te.getFrame();
        if (frame == null) {
            frame = FrameManager.getDefaultFrame();  // Clicked over name panel, not a specific frame
        }
        final Range range = frame.getCurrentRange();
        final String chrom = range.getChr();
        final int chromStart = (int) frame.getChromosomePosition(me);

        // Change track height by attribute
        JMenu groupMenu = new JMenu("Group alignments by");
        ButtonGroup group = new ButtonGroup();

        AlignmentTrack.GroupOption[] groupOptions = {
                AlignmentTrack.GroupOption.NONE, AlignmentTrack.GroupOption.STRAND,
                AlignmentTrack.GroupOption.FIRST_OF_PAIR_STRAND, AlignmentTrack.GroupOption.SAMPLE,
                AlignmentTrack.GroupOption.LIBRARY, AlignmentTrack.GroupOption.READ_GROUP,
                AlignmentTrack.GroupOption.MATE_CHROMOSOME,
                AlignmentTrack.GroupOption.PAIR_ORIENTATION, AlignmentTrack.GroupOption.CHIMERIC,
                AlignmentTrack.GroupOption.SUPPLEMENTARY, AlignmentTrack.GroupOption.REFERENCE_CONCORDANCE,
                AlignmentTrack.GroupOption.MOVIE, AlignmentTrack.GroupOption.ZMW, AlignmentTrack.GroupOption.READ_ORDER,
                AlignmentTrack.GroupOption.LINKED, AlignmentTrack.GroupOption.PHASE,
                AlignmentTrack.GroupOption.MAPPING_QUALITY,
                AlignmentTrack.GroupOption.SELECTED,
                AlignmentTrack.GroupOption.DUPLICATE
        };

        for (final AlignmentTrack.GroupOption option : groupOptions) {
            JCheckBoxMenuItem mi = new JCheckBoxMenuItem(option.label);
            mi.setSelected(renderOptions.getGroupByOption() == option);
            mi.addActionListener(aEvt -> groupAlignments(option, null, null));
            groupMenu.add(mi);
            group.add(mi);
        }

        JCheckBoxMenuItem tagOption = new JCheckBoxMenuItem("tag");
        tagOption.addActionListener(aEvt -> {
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getGroupByTag());
            if (tag != null) {
                if (tag.trim().length() > 0) {
                    groupAlignments(AlignmentTrack.GroupOption.TAG, tag, null);
                } else {
                    groupAlignments(AlignmentTrack.GroupOption.NONE, null, null);
                }
            }

        });
        tagOption.setSelected(renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.TAG);
        groupMenu.add(tagOption);
        group.add(tagOption);

//        Range oldGroupByPos = renderOptions.getGroupByPos();
//        if (oldGroupByPos != null && renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.BASE_AT_POS) { // already sorted by the base at a position
//            JCheckBoxMenuItem oldGroupByPosOption = new JCheckBoxMenuItem("base at " + oldGroupByPos.getChr() +
//                    ":" + Globals.DECIMAL_FORMAT.format(1 + oldGroupByPos.getStart()));
//            groupMenu.add(oldGroupByPosOption);
//            oldGroupByPosOption.setSelected(true);
//        }
//
//        if (renderOptions.getGroupByOption() != AlignmentTrack.GroupOption.BASE_AT_POS || oldGroupByPos == null ||
//                !oldGroupByPos.getChr().equals(chrom) || oldGroupByPos.getStart() != chromStart) { // not already sorted by this position
        JCheckBoxMenuItem newGroupByPosOption = new JCheckBoxMenuItem("base at " + chrom +
                ":" + Globals.DECIMAL_FORMAT.format(1 + chromStart));
        newGroupByPosOption.addActionListener(aEvt -> {
            Range groupByPos = new Range(chrom, chromStart, chromStart + 1);
            groupAlignments(AlignmentTrack.GroupOption.BASE_AT_POS, null, groupByPos);
        });
        groupMenu.add(newGroupByPosOption);
        group.add(newGroupByPosOption);
        // }

        JCheckBoxMenuItem newGroupByInsOption = new JCheckBoxMenuItem("insertion at " + chrom +
                ":" + Globals.DECIMAL_FORMAT.format(1 + chromStart));
        newGroupByInsOption.addActionListener(aEvt -> {
            Range groupByPos = new Range(chrom, chromStart, chromStart + 1);
            groupAlignments(AlignmentTrack.GroupOption.INSERTION_AT_POS, null, groupByPos);
        });
        groupMenu.add(newGroupByInsOption);
        group.add(newGroupByInsOption);


        groupMenu.add(new Separator());
        JCheckBoxMenuItem invertGroupNameSortingOption = new JCheckBoxMenuItem("Reverse group order");
        invertGroupNameSortingOption.setSelected(renderOptions.isInvertGroupSorting());
        invertGroupNameSortingOption.addActionListener(aEvt -> {
            renderOptions.setInvertGroupSorting(!renderOptions.isInvertGroupSorting());
            alignmentTrack.packAlignments();
            alignmentTrack.repaint();
        });
        groupMenu.add(invertGroupNameSortingOption);

        JCheckBoxMenuItem groupAllOption = new JCheckBoxMenuItem("Group all tracks");
        groupAllOption.setSelected(alignmentTrack.getPreferences().getAsBoolean(SAM_GROUP_ALL));
        groupAllOption.addActionListener(aEvt -> {
            alignmentTrack.getPreferences().put(SAM_GROUP_ALL, groupAllOption.getState());
        });
        groupMenu.add(groupAllOption);

        add(groupMenu);
    }

    private void groupAlignments(AlignmentTrack.GroupOption option, String tag, Range pos) {

        if (alignmentTrack.getPreferences().getAsBoolean(SAM_GROUP_ALL)) {
            for (AlignmentTrack t : IGV.getInstance().getAlignmentTracks()) {
                t.groupAlignments(option, tag, pos);
            }
        } else {
            alignmentTrack.groupAlignments(option, tag, pos);
        }
    }


    /**
     * Sort menu
     */
    void addSortMenuItem() {

        JMenu sortMenu = new JMenu("Sort alignments by");
        //LinkedHashMap is supposed to preserve order of insertion for iteration
        Map<String, SortOption> mappings = new LinkedHashMap<>();

        mappings.put("start location", SortOption.START);
        mappings.put("read strand", SortOption.STRAND);
        mappings.put("first-of-pair strand", SortOption.FIRST_OF_PAIR_STRAND);
        mappings.put("base", SortOption.NUCLEOTIDE);
        mappings.put("mapping quality", SortOption.QUALITY);
        mappings.put("sample", SortOption.SAMPLE);
        mappings.put("read group", SortOption.READ_GROUP);
        mappings.put("read order", SortOption.READ_ORDER);
        mappings.put("read name", SortOption.READ_NAME);
        mappings.put("aligned read length", SortOption.ALIGNED_READ_LENGTH);
        mappings.put("left clip", SortOption.LEFT_CLIP);
        mappings.put("right clip", SortOption.RIGHT_CLIP);

        if (dataManager.isPairedEnd()) {
            mappings.put("insert size", SortOption.INSERT_SIZE);
            mappings.put("chromosome of mate", SortOption.MATE_CHR);
        }

        for (Map.Entry<String, SortOption> el : mappings.entrySet()) {
            JMenuItem mi = new JMenuItem(el.getKey());
            mi.addActionListener(aEvt -> {
                final SortOption option = el.getValue();
                renderOptions.setSortOption(option);
                alignmentTrack.sortAlignmentTracks(option, null, renderOptions.isInvertSorting());
            });
            sortMenu.add(mi);
        }

        JMenuItem tagOption = new JMenuItem("tag");
        tagOption.addActionListener(aEvt -> {
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getSortByTag());
            if (tag != null && tag.trim().length() > 0) {
                renderOptions.setSortByTag(tag);
                renderOptions.setSortOption((SortOption.TAG));
                alignmentTrack.sortAlignmentTracks(SortOption.TAG, tag, renderOptions.isInvertSorting());
            }
        });
        sortMenu.add(tagOption);

        sortMenu.add(new Separator());
        JCheckBoxMenuItem invertGroupNameSortingOption = new JCheckBoxMenuItem("reverse sorting");
        invertGroupNameSortingOption.setSelected(renderOptions.isInvertSorting());
        invertGroupNameSortingOption.addActionListener(aEvt -> {
            final boolean updatedInvertSorting = !renderOptions.isInvertSorting();
            renderOptions.setInvertSorting(updatedInvertSorting);
            alignmentTrack.sortAlignmentTracks(renderOptions.getSortOption(), renderOptions.getSortByTag(), updatedInvertSorting);
        });
        sortMenu.add(invertGroupNameSortingOption);
        add(sortMenu);
    }

    public void addFilterMenuItem() {
        JMenu filterMenu = new JMenu("Filter alignments by");
        JMenuItem mi = new JMenuItem("mapping quality");
        mi.addActionListener(aEvt -> {
            // TODO -- use current value for default
            String defString = PreferencesManager.getPreferences().get(SAM_QUALITY_THRESHOLD);
            if (defString == null) defString = "";
            String mqString = MessageUtils.showInputDialog("Minimum mapping quality: ", defString);
            try {
                int mq = Integer.parseInt(mqString);
                // TODO do something with this
                //System.out.println(mq);
            } catch (NumberFormatException e) {
                MessageUtils.showMessage("Mapping quality must be an integer");
            }
        });
        filterMenu.add(mi);
        add(filterMenu);
    }

    private JRadioButtonMenuItem getColorMenuItem(String label, final AlignmentTrack.ColorOption option) {
        JRadioButtonMenuItem mi = new JRadioButtonMenuItem(label);
        mi.setSelected(renderOptions.getColorOption() == option);
        mi.addActionListener(aEvt -> {
            alignmentTrack.setColorOption(option);
            alignmentTrack.repaint();
        });
        return mi;
    }

    private JRadioButtonMenuItem getBasemodColorMenuItem(String label, final AlignmentTrack.ColorOption option, boolean groupByStrand, String filter) {
        JRadioButtonMenuItem mi = new JRadioButtonMenuItem(label);
        mi.setSelected(renderOptions.getColorOption() == option);
        mi.addActionListener(aEvt -> {
            alignmentTrack.setColorOption(option);
            renderOptions.setBasemodFilter(filter == null ? null : new BaseModficationFilter(filter));
            if (groupByStrand) {
                alignmentTrack.groupAlignments(AlignmentTrack.GroupOption.FIRST_OF_PAIR_STRAND, null, null);
            } else {
                alignmentTrack.repaint();
            }
        });

        return mi;
    }

    void addColorByMenuItem() {
        // Change track height by attribute
        JMenu colorMenu = new JMenu("Color alignments by");

        ButtonGroup group = new ButtonGroup();

        Map<String, AlignmentTrack.ColorOption> mappings = new LinkedHashMap<>();

        mappings.put("none", AlignmentTrack.ColorOption.NONE);

        if (dataManager.hasYCTags()) {
            mappings.put("YC tag", AlignmentTrack.ColorOption.YC_TAG);
        }

        if (dataManager.isPairedEnd()) {
            mappings.put("insert size", AlignmentTrack.ColorOption.INSERT_SIZE);
            mappings.put("pair orientation", AlignmentTrack.ColorOption.PAIR_ORIENTATION);
            mappings.put("insert size and pair orientation", AlignmentTrack.ColorOption.UNEXPECTED_PAIR);
        }

        mappings.put("read strand", AlignmentTrack.ColorOption.READ_STRAND);

        if (dataManager.isPairedEnd()) {
            mappings.put("first-of-pair strand", AlignmentTrack.ColorOption.FIRST_OF_PAIR_STRAND);
        }

        mappings.put("read group", AlignmentTrack.ColorOption.READ_GROUP);

        if (dataManager.isPairedEnd()) {
            mappings.put("read order", AlignmentTrack.ColorOption.READ_ORDER);
        }

        mappings.put("sample", AlignmentTrack.ColorOption.SAMPLE);
        mappings.put("library", AlignmentTrack.ColorOption.LIBRARY);
        mappings.put("movie", AlignmentTrack.ColorOption.MOVIE);
        mappings.put("ZMW", AlignmentTrack.ColorOption.ZMW);
        mappings.put("split", AlignmentTrack.ColorOption.SPLIT);

        for (Map.Entry<String, AlignmentTrack.ColorOption> el : mappings.entrySet()) {
            JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
            colorMenu.add(mi);
            group.add(mi);
        }

        JRadioButtonMenuItem tagOption = new JRadioButtonMenuItem("tag");
        tagOption.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.TAG);
        tagOption.addActionListener(aEvt -> {
            alignmentTrack.setColorOption(AlignmentTrack.ColorOption.TAG);
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getColorByTag());
            if (tag != null && tag.trim().length() > 0) {
                alignmentTrack.setColorByTag(tag);
                alignmentTrack.repaint();
            }
        });
        colorMenu.add(tagOption);
        group.add(tagOption);

        colorMenu.add(getBisulfiteContextMenuItem(group));

        // Base modifications
        JRadioButtonMenuItem bmMenuItem;

        Set<String> allModifications = dataManager.getAllBaseModificationKeys().stream().map(BaseModificationKey::getModification).collect(Collectors.toSet());
        final int modificationCount = allModifications.size();
        if (modificationCount > 0) {
            BaseModficationFilter filter = renderOptions.getBasemodFilter();
            boolean groupByStrand = alignmentTrack.getPreferences().getAsBoolean(BASEMOD_GROUP_BY_STRAND);
            colorMenu.addSeparator();
            if (modificationCount > 1) {
                bmMenuItem = getBasemodColorMenuItem("base modification (all)", AlignmentTrack.ColorOption.BASE_MODIFICATION, groupByStrand, null);
                bmMenuItem.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.BASE_MODIFICATION && filter == null);
                colorMenu.add(bmMenuItem);
                group.add(bmMenuItem);
            }
            for (String m : allModifications) {
                String name = BaseModificationUtils.modificationName(m);
                bmMenuItem = getBasemodColorMenuItem("base modification (" + name + ")", AlignmentTrack.ColorOption.BASE_MODIFICATION, groupByStrand, m);
                bmMenuItem.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.BASE_MODIFICATION && (filter != null && filter.pass(m)));
                colorMenu.add(bmMenuItem);
                group.add(bmMenuItem);
            }


            colorMenu.addSeparator();
            if (modificationCount > 1) {
                bmMenuItem = getBasemodColorMenuItem("base modification 2-color (all)", AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR, groupByStrand, null);
                bmMenuItem.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR && filter == null);
                colorMenu.add(bmMenuItem);
                group.add(bmMenuItem);
            }
            for (String m : allModifications) {
                String name = BaseModificationUtils.modificationName(m);
                bmMenuItem = getBasemodColorMenuItem("base modification 2-color (" + name + ")", AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR, groupByStrand, m);
                bmMenuItem.setSelected(renderOptions.getColorOption() ==
                        AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR &&
                        (filter != null && filter.pass(m)));
                colorMenu.add(bmMenuItem);
                group.add(bmMenuItem);
            }
        }


        // SMRT kinetics
        if (alignmentTrack.getPreferences().getAsBoolean(SMRT_KINETICS_SHOW_OPTIONS)) {
            // Show additional options to help visualize SMRT kinetics data
            mappings.clear();
            mappings.put("SMRT subread IPD", AlignmentTrack.ColorOption.SMRT_SUBREAD_IPD);
            mappings.put("SMRT subread PW", AlignmentTrack.ColorOption.SMRT_SUBREAD_PW);
            mappings.put("SMRT CCS fwd-strand aligned IPD", AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD);
            mappings.put("SMRT CCS fwd-strand aligned PW", AlignmentTrack.ColorOption.SMRT_CCS_FWD_PW);
            mappings.put("SMRT CCS rev-strand aligned IPD", AlignmentTrack.ColorOption.SMRT_CCS_REV_IPD);
            mappings.put("SMRT CCS rev-strand aligned PW", AlignmentTrack.ColorOption.SMRT_CCS_REV_PW);
            colorMenu.addSeparator();
            for (Map.Entry<String, AlignmentTrack.ColorOption> el : mappings.entrySet()) {
                JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
                colorMenu.add(mi);
                group.add(mi);
            }
        }

        add(colorMenu);

    }

    void addShadeAlignmentsMenuItem() {
        JMenu shadeMenu = new JMenu("Shade alignments by");

        for (AlignmentTrack.ShadeAlignmentsOption option : AlignmentTrack.ShadeAlignmentsOption.values()) {
            JRadioButtonMenuItem mi = new JRadioButtonMenuItem(option.label);
            mi.setSelected(renderOptions.getShadeAlignmentsOption() == option);
            mi.addActionListener(aEvt -> {
                alignmentTrack.setShadeAlignmentsOptions(option);
                alignmentTrack.repaint();
            });

            shadeMenu.add(mi);
        }
        add(shadeMenu);
    }

    void addPackMenuItem() {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Re-pack alignments");
        item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            IGV.getInstance().packAlignmentTracks();
            alignmentTrack.repaint();
        }));

        add(item);
    }

    void addCopyToClipboardItem(final TrackClickEvent te, Alignment alignment) {

        final MouseEvent me = te.getMouseEvent();
        JMenuItem item = new JMenuItem("Copy read details");
        final ReferenceFrame frame = te.getFrame();
        if (frame == null) {
            item.setEnabled(false);
        } else {
            final double location = frame.getChromosomePosition(me);

            // Change track height by attribute
            item.addActionListener(aEvt -> copyToClipboard(te, alignment, location, me.getX()));
            if (alignment == null) {
                item.setEnabled(false);
            }
        }
        add(item);
    }


    void addViewAsPairsMenuItem() {
        final JMenuItem item = new JCheckBoxMenuItem("View as pairs");
        item.setSelected(renderOptions.isViewPairs());
        item.addActionListener(aEvt -> {
            boolean viewAsPairs = item.isSelected();
            alignmentTrack.setViewAsPairs(viewAsPairs);
        });
        item.setEnabled(dataManager.isPairedEnd());
        add(item);
    }

    void addGoToMate(final TrackClickEvent te, Alignment clickedAlignment) {
        JMenuItem item = new JMenuItem("Go to mate");
        addActionIfMatesAreMapped(te, clickedAlignment, item, this::gotoMate);
        add(item);
    }

    void showMateRegion(final TrackClickEvent te, Alignment clickedAlignment) {
        JMenuItem item = new JMenuItem("View mate region in split screen");
        addActionIfMatesAreMapped(te, clickedAlignment, item, this::splitScreenMate);
        add(item);
    }

    private void addActionIfMatesAreMapped(final TrackClickEvent te, final Alignment clickedAlignment, final JMenuItem item,
                                           final BiConsumer<ReferenceFrame, Alignment> action) {
        final ReferenceFrame frame = te.getFrame();
        if (frame == null) {
            item.setEnabled(false);
        } else {
            final Alignment alignment = clickedAlignment.getSpecificAlignment(te.getChromosomePosition()); // handle PairedAlignments and the like
            item.addActionListener(aEvt -> action.accept(frame, alignment));
            if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                item.setEnabled(false);
            }
        }
    }

    void addClearSelectionsMenuItem() {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Clear selections");
        item.addActionListener(aEvt -> {
            alignmentTrack.getSelectedReadNames().clear();
            alignmentTrack.repaint();
        });
        add(item);
    }

    JMenuItem addShowAllBasesMenuItem() {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show all bases");

        if (renderOptions.getColorOption() == AlignmentTrack.ColorOption.BISULFITE || renderOptions.getColorOption() == AlignmentTrack.ColorOption.NOMESEQ) {
            //    item.setEnabled(false);
        } else {
            item.setSelected(renderOptions.isShowAllBases());
        }
        item.addActionListener(aEvt -> {
            renderOptions.setShowAllBases(item.isSelected());
            alignmentTrack.repaint();
        });
        add(item);
        return item;
    }

    JMenuItem addShowMismatchesMenuItem() {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show mismatched bases");
        item.setSelected(renderOptions.isShowMismatches());
        item.addActionListener(aEvt -> {
            renderOptions.setShowMismatches(item.isSelected());
            alignmentTrack.repaint();
        });
        add(item);
        return item;
    }


    void addInsertSizeMenuItem() {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Set insert size options ...");
        item.addActionListener(aEvt -> {
            InsertSizeSettingsDialog dlg = new InsertSizeSettingsDialog(IGV.getInstance().getMainFrame(), renderOptions);
            dlg.setModal(true);
            dlg.setVisible(true);
            if (!dlg.isCanceled()) {
                renderOptions.setComputeIsizes(dlg.isComputeIsize());
                renderOptions.setMinInsertSizePercentile(dlg.getMinPercentile());
                renderOptions.setMaxInsertSizePercentile(dlg.getMaxPercentile());
                if (renderOptions.computeIsizes) {
                    dataManager.updatePEStats(renderOptions);
                }

                renderOptions.setMinInsertSize(dlg.getMinThreshold());
                renderOptions.setMaxInsertSize(dlg.getMaxThreshold());
                alignmentTrack.repaint();
            }
        });


        item.setEnabled(dataManager.isPairedEnd());
        add(item);
    }

    void addShadeBaseByMenuItem() {

        final JMenuItem item = new JCheckBoxMenuItem("Shade base by quality");
        item.setSelected(renderOptions.getShadeBasesOption());
        item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            renderOptions.setShadeBasesOption(item.isSelected());
            alignmentTrack.repaint();
        }));
        add(item);
    }

    void addShowItems() {

        final CoverageTrack coverageTrack = alignmentTrack.getCoverageTrack();
        if (coverageTrack != null) {
            final JMenuItem item = new JCheckBoxMenuItem("Show Coverage Track");
            item.setSelected(coverageTrack.isVisible());
            item.setEnabled(!coverageTrack.isRemoved());
            item.addActionListener(aEvt -> {
                coverageTrack.setVisible(item.isSelected());
                IGV.getInstance().repaint(Arrays.asList(coverageTrack));

            });
            add(item);
        }

        final SpliceJunctionTrack spliceJunctionTrack = alignmentTrack.getSpliceJunctionTrack();
        if (spliceJunctionTrack != null) {
            final JMenuItem item = new JCheckBoxMenuItem("Show Splice Junction Track");
            item.setSelected(spliceJunctionTrack.isVisible());
            item.setEnabled(!spliceJunctionTrack.isRemoved());
            item.addActionListener(aEvt -> {
                alignmentTrack.setVisible(item.isSelected());
                IGV.getInstance().repaint(Arrays.asList(spliceJunctionTrack));

            });
            add(item);
        }

        final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
        alignmentItem.setSelected(true);
        alignmentItem.addActionListener(e -> {
            alignmentTrack.setVisible(alignmentItem.isSelected());
            IGV.getInstance().repaint(Arrays.asList(alignmentTrack));
        });
        // Disable if this is the only visible track
        if (!((coverageTrack != null && coverageTrack.isVisible()) ||
                (spliceJunctionTrack != null && spliceJunctionTrack.isVisible()))) {
            alignmentItem.setEnabled(false);
        }

        add(alignmentItem);
    }


    void addCopySequenceItems(final TrackClickEvent te) {

        final JMenuItem item = new JMenuItem("Copy read sequence");
        add(item);
        final Alignment alignment = getSpecificAlignment(te);
        if (alignment == null) {
            item.setEnabled(false);
            return;
        }
        final String seq = alignment.getReadSequence();
        if (seq == null) {
            item.setEnabled(false);
            return;
        }
        item.addActionListener(aEvt -> StringUtils.copyTextToClipboard(seq));

        /* Add a "Copy left clipped sequence" item if there is  left clipping. */
        int minimumBlatLength = BlatClient.MINIMUM_BLAT_LENGTH;
        ClippingCounts clipping = alignment.getClippingCounts();
        if (clipping.getLeftSoft() > 0) {
            String lcSeq = getClippedSequence(alignment.getReadSequence(), 0, clipping.getLeftSoft());
            final JMenuItem lccItem = new JMenuItem("Copy left-clipped sequence");
            add(lccItem);
            lccItem.addActionListener(aEvt -> StringUtils.copyTextToClipboard(lcSeq));
        }

        /* Add a "Copy right clipped sequence" item if there is  right clipping. */
        if (clipping.getRightSoft() > 0) {
            int seqLength = seq.length();
            String rcSeq = getClippedSequence(
                    alignment.getReadSequence(),
                    seqLength - clipping.getRightSoft(),
                    seqLength);

            final JMenuItem rccItem = new JMenuItem("Copy right-clipped sequence");
            add(rccItem);
            rccItem.addActionListener(aEvt -> StringUtils.copyTextToClipboard(rcSeq));
        }
    }


    void addBlatItem(final TrackClickEvent te) {
        // Change track height by attribute
        final JMenuItem item = new JMenuItem("BLAT read sequence");
        add(item);

        final Alignment alignment = getSpecificAlignment(te);
        if (alignment == null) {
            item.setEnabled(false);
            return;
        }

        final String seq = alignment.getReadSequence();
        if (seq == null || seq.equals("*")) {
            item.setEnabled(false);
            return;

        }

        item.addActionListener(aEvt -> {
            String blatSeq = alignment.getReadStrand() == Strand.NEGATIVE ?
                    SequenceTrack.getReverseComplement(seq) : seq;
            BlatClient.doBlatQuery(blatSeq, alignment.getReadName());
        });

    }

    void addBlatClippingItems(final TrackClickEvent te) {
        final Alignment alignment = getSpecificAlignment(te);
        if (alignment == null) {
            return;
        }

        int minimumBlatLength = BlatClient.MINIMUM_BLAT_LENGTH;
        ClippingCounts clipping = alignment.getClippingCounts();

        /* Add a "BLAT left clipped sequence" item if there is significant left clipping. */
        if (clipping.getLeftSoft() > minimumBlatLength) {
            String lcSeq = getClippedSequence(alignment.getReadSequence(), 0, clipping.getLeftSoft());
            String blatSeq = alignment.isNegativeStrand() ?
                    SequenceTrack.getReverseComplement(lcSeq) :
                    lcSeq;
            final JMenuItem lcbItem = new JMenuItem("BLAT left-clipped sequence");
            add(lcbItem);
            lcbItem.addActionListener(aEvt ->
                    BlatClient.doBlatQuery(blatSeq, alignment.getReadName() + " - left clip")
            );
        }
        /* Add a "BLAT right clipped sequence" item if there is significant right clipping. */
        if (clipping.getRightSoft() > minimumBlatLength) {

            String seq = alignment.getReadSequence();
            int seqLength = seq.length();
            String rcSeq = getClippedSequence(
                    alignment.getReadSequence(),
                    seqLength - clipping.getRightSoft(),
                    seqLength);
            String blatSeq = alignment.isNegativeStrand() ?
                    SequenceTrack.getReverseComplement(rcSeq) :
                    rcSeq;

            final JMenuItem rcbItem = new JMenuItem("BLAT right-clipped sequence");
            add(rcbItem);
            rcbItem.addActionListener(aEvt ->
                    BlatClient.doBlatQuery(blatSeq, alignment.getReadName() + " - right clip")
            );

        }
    }

    private String getClippedSequence(String readSequence, int i, int i2) {
        if (readSequence == null || readSequence.equals("*")) {
            return "*";
        }
        return readSequence.substring(i, i2);
    }

    void addExtViewItem(final TrackClickEvent te) {
        // Change track height by attribute
        final JMenuItem item = new JMenuItem("ExtView");
        add(item);

        final Alignment alignment = alignmentTrack.getAlignmentAt(te);
        if (alignment == null) {
            item.setEnabled(false);
            return;
        }

        final String seq = alignment.getReadSequence();
        if (seq == null) {
            item.setEnabled(false);
            return;

        }

        item.addActionListener(aEvt -> ExtendViewClient.postExtendView(alignment));

    }

    /**
     * Add all menu items that link alignments by tag or readname.  These are mutually exclusive.  The
     * list includes 2 items for 10X "Loupe link-read" style views, a supplementary alignment option,
     * and linking by arbitrary tag.
     */
    void addLinkedReadItems() {

        final JCheckBoxMenuItem supplementalItem = new JCheckBoxMenuItem("Link supplementary alignments");
        supplementalItem.setSelected(alignmentTrack.isLinkedReads() && "READNAME".equals(renderOptions.getLinkByTag()));
        supplementalItem.addActionListener(aEvt -> {
            boolean linkedReads = supplementalItem.isSelected();
            setLinkByTag(linkedReads, "READNAME");
        });
        add(supplementalItem);

        String linkedTagsString = PreferencesManager.getPreferences().get(SAM_LINK_BY_TAGS);
        if (linkedTagsString != null) {
            String[] t = Globals.commaPattern.split(linkedTagsString);
            for (String tag : t) {
                if (tag.length() > 0) {
                    add(linkedReadItem(tag));
                }
            }
        }

        final JMenuItem linkByTagItem = new JMenuItem("Link by tag...");
        linkByTagItem.addActionListener(aEvt -> {
            String tag = MessageUtils.showInputDialog("Link by tag:");
            if (tag != null) {
                setLinkByTag(true, tag);
                String linkedTags = PreferencesManager.getPreferences().get(SAM_LINK_BY_TAGS);
                if (linkedTags == null) {
                    linkedTags = tag;
                } else {
                    linkedTags += "," + tag;
                }
                PreferencesManager.getPreferences().put(SAM_LINK_BY_TAGS, linkedTags);
            }
        });
        add(linkByTagItem);
    }

    private JCheckBoxMenuItem linkedReadItem(String tag) {
        final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Link by " + tag);
        item.setSelected(!alignmentTrack.isLinkedReadView() && alignmentTrack.isLinkedReads() && tag.equals(renderOptions.getLinkByTag()));
        item.addActionListener(aEvt -> {
            boolean linkedReads = item.isSelected();
            if ("BX".equals(tag) || "MI".equals(tag)) {
                alignmentTrack.setLinkedReadView(linkedReads, tag);
            } else {
                setLinkByTag(linkedReads, tag);
            }
        });
        return item;
    }


    private void addInsertionItems(AlignmentBlock insertion) {

        final JMenuItem item = new JMenuItem("Copy insert sequence");
        add(item);
        item.addActionListener(aEvt -> StringUtils.copyTextToClipboard(insertion.getBases().getString()));

        if (insertion.getBases() != null && insertion.getBases().length >= 10) {
            final JMenuItem blatItem = new JMenuItem("BLAT insert sequence");
            add(blatItem);
            blatItem.addActionListener(aEvt -> {
                String blatSeq = insertion.getBases().getString();
                BlatClient.doBlatQuery(blatSeq, "BLAT insert sequence");
            });
        }
    }

    void addThirdGenItems(Alignment clickedAlignment, final TrackClickEvent tce) {

        // Linked read items -- mostly for 3rd gen but might also be relevant to 10X and other linked read assays
        addLinkedReadItems();

        //Supplementary/chimeric items, only if the read has an SA tag;
        addShowChimericRegions(alignmentTrack, tce, clickedAlignment);
        addShowDiagram(tce, clickedAlignment);

    }

    void addSBXItems(Alignment clickedAlignment, final TrackClickEvent tce) {

        addSeparator();

        final JMenuItem item = new JCheckBoxMenuItem("INDEL coloring uses grey (SBX)");
        item.setSelected(renderOptions.isIndelQualSbx());
        item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            renderOptions.setIndelQualSbx(item.isSelected());
            alignmentTrack.repaint();
        }));
        add(item);

        final JMenuItem item2 = new JCheckBoxMenuItem("Simplex tail coloring (SBX)");
        item2.setSelected(renderOptions.isTailQualSbx());
        item2.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            renderOptions.setTailQualSbx(item2.isSelected());
            alignmentTrack.repaint();
        }));
        add(item2);

        final JMenuItem item3 = new JCheckBoxMenuItem("Hide simplex tails (SBX)");
        item3.setSelected(renderOptions.isHideTailSbx());
        item3.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            renderOptions.setHideTailSbx(item3.isSelected());
            alignmentTrack.repaint();
        }));
        add(item3);
    }

    /**
     * Copy the contents of the popup text to the system clipboard.
     */
    private void copyToClipboard(final TrackClickEvent e, Alignment alignment, double location, int mouseX) {

        if (alignment != null) {
            final String clipboardString = alignment.getClipboardString(location, mouseX)
                    .replace("<b>", "")
                    .replace("</b>", "")
                    .replace("<br>", "\n")
                    .replace("<br/>", "\n")
                    .replace("<hr>", "\n------------------\n")
                    .replace("<hr/>", "\n------------------\n");
            StringSelection stringSelection = new StringSelection(clipboardString);
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            clipboard.setContents(stringSelection, null);
        }
    }


    /**
     * Jump to the mate region
     */
    private void gotoMate(final ReferenceFrame frame, Alignment alignment) {
        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                alignmentTrack.setSelectedAlignment(alignment);

                String chr = mate.getChr();
                int start = mate.start - 1;

                // Don't change scale
                double range = frame.getEnd() - frame.getOrigin();
                int newStart = (int) Math.max(0, (start + (alignment.getEnd() - alignment.getStart()) / 2 - range / 2));
                int newEnd = newStart + (int) range;
                frame.jumpTo(chr, newStart, newEnd);
                AlignmentTrack.sortSelectedReadsToTheTop(alignmentTrack.getSelectedReadNames().keySet());
                frame.recordHistory();
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }


    /**
     * Split the screen so the current view and mate region are side by side.
     * Need a better name for this method.
     */
    private void splitScreenMate(ReferenceFrame frame, Alignment alignment) {
        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {
                alignmentTrack.setSelectedAlignment(alignment);
                FrameManager.addNewLociToFrames(frame, List.of(mate), alignmentTrack.getSelectedReadNames().keySet());
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }


    /**
     * Get the most "specific" alignment at the specified location.  Specificity refers to the smallest alignment
     * in a group that contains the location (i.e. if a group of linked alignments overlap take the smallest one).
     *
     * @param te
     * @return
     */
    Alignment getSpecificAlignment(TrackClickEvent te) {
        Alignment alignment = alignmentTrack.getAlignmentAt(te);
        if (alignment != null) {
            alignment = alignment.getSpecificAlignment(te.getChromosomePosition());
        }
        return alignment;
    }


    /**
     * Link alignments by arbitrary tag, without the extra settings applied to link-read-view
     *
     * @param linkReads
     * @param tag
     */
    private void setLinkByTag(boolean linkReads, String tag) {
        if (alignmentTrack.isLinkedReadView()) {
            alignmentTrack.undoLinkedReadView();
        }
        if (linkReads) {
            renderOptions.setLinkByTag(tag);
            if (renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.NONE) {
                renderOptions.setGroupByOption(AlignmentTrack.GroupOption.LINKED);
            }
        } else {
            renderOptions.setLinkByTag(null);
            if (renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.LINKED) {
                renderOptions.setGroupByOption(AlignmentTrack.GroupOption.NONE);
            }
        }
        renderOptions.setLinkedReads(linkReads);
        alignmentTrack.packAlignments();
        alignmentTrack.repaint();
    }


    /**
     * Listener for deselecting one component when another is selected
     */
    private static class Deselector implements ActionListener {

        private final JMenuItem toDeselect;
        private final JMenuItem parent;

        Deselector(JMenuItem parent, JMenuItem toDeselect) {
            this.parent = parent;
            this.toDeselect = toDeselect;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            if (this.parent.isSelected()) {
                this.toDeselect.setSelected(false);
            }
        }
    }

}

