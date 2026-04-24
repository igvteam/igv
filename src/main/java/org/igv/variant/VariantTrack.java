//chr2:128,565,093-128,565,156

package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.tribble.TribbleException;
import org.igv.Globals;
import org.igv.event.IGVEventObserver;
import org.igv.feature.FeatureUtils;
import org.igv.feature.PackedFeature;
import org.igv.jbrowse.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.GraphicUtils;
import org.igv.sample.SampleGroup;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.StringUtils;
import org.igv.variant.vcf.MateVariant;
import org.igv.variant.vcf.VCFVariant;
import org.w3c.dom.Element;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static org.igv.prefs.Constants.DEFAULT_VISIBILITY_WINDOW;
import static org.igv.prefs.Constants.VARIANT_COLOR_BY_ALLELE_FREQ;

/**
 * @author Jesse Whitworth, Jim Robinson, Fabien Campagne
 */

public class VariantTrack extends FeatureTrack implements IGVEventObserver {


    public static final int DEFAULT_MAX_HEIGHT = 300;
    private static Logger log = LogManager.getLogger(VariantTrack.class);

    static final DecimalFormat numFormat = new DecimalFormat("#.###");

    private static final Color CIRC_VIEW_DEFAULT_COLOR = new Color(27, 192, 249);
    private static final int GROUP_BORDER_WIDTH = 0;
    private static final Color BAND1_COLOR = new Color(245, 245, 245);
    private static final Color BAND2_COLOR = Globals.isDarkMode() ? new Color(200, 200, 200) : Color.white;

    private final static int DEFAULT_EXPANDED_GENOTYPE_HEIGHT = 15;
    private final static int DEFAULT_EXPANDED_VARIANT_HEIGHT = 25;
    private final static int DEFAULT_SQUISHED_VARIANT_HEIGHT = 6;
    private final static int MAX_FILTER_LINES = 15;
    private final static int WG_TRACK_HEIGHT = 40;
    private final int DEFAULT_SQUISHED_HEIGHT = 4;


    // TODO -- this needs to be settable
    public static int METHYLATION_MIN_BASE_COUNT = 10;
    private transient Rectangle lastClipBounds;

    public static boolean isVCF(String format) {
        return (format.equals("vcf3") ||
                format.equals("vcf4") ||
                format.equals("vcf") ||
                format.equals("bcf") ||
                format.equals("gvcf"));
    }


    private VariantRenderer renderer;

    /**
     * When this flag is true, we have detected that the VCF file contains the FORMAT MR column representing
     * methylation data. This will enable the "Color By/Methylation Rate" menu item.
     */
    private boolean enableMethylationRateSupport;

    /**
     * Top (y) position of this track.  This is updated whenever the track is drawn.
     */

    private boolean showGenotypes = true;

    /**
     * The height of a single row in in squished mode
     */
    private int squishedHeight = DEFAULT_SQUISHED_HEIGHT;

    /**
     * Current coloring option
     */
    private ColorMode genotypeColorMode = ColorMode.GENOTYPE;


    private ColorMode siteColorMode;

    /**
     * When true, variants that are marked filtering are not drawn.
     */
    private boolean hideFiltered = false;

    /**
     * The currently selected variant.  This is a transient variable, set only while the popup menu is up.
     */
    private Variant selectedVariant;

    /**
     * Transient list to keep track of the vertical bounds of each sample.  Set when rendering names, used to
     * select correct sample for popup text.  We use a list and linear lookup for now, some sort of tree structure
     * would be faster.
     */
    private List<SampleBounds> sampleBounds = new ArrayList<SampleBounds>();

    public void setRenderer(VariantRenderer renderer) {
        this.renderer = renderer;
    }

    public VariantTrack() {
    }

    public VariantTrack(String name, FeatureSource source) {
        this(null, source, Collections.<String>emptyList(), false);
        this.setName(name);
    }

    public VariantTrack(ResourceLocator locator, FeatureSource source, List<String> samples,
                        boolean enableMethylationRateSupport) {
        super(locator, source);

        IGVPreferences prefMgr = PreferencesManager.getPreferences();

        this.renderer = new VariantRenderer(this);

        this.enableMethylationRateSupport = enableMethylationRateSupport;
        if (enableMethylationRateSupport) {
            // also set the default color mode to Methylation rate:
            genotypeColorMode = ColorMode.METHYLATION_RATE;
        }

        this.siteColorMode = prefMgr.getAsBoolean(VARIANT_COLOR_BY_ALLELE_FREQ) ?
                ColorMode.ALLELE_FREQUENCY :
                ColorMode.ALLELE_FRACTION;

        this.initSamples(samples);

        setDisplayMode(DisplayMode.EXPANDED);

        int sampleCount = samples.size();
        final int groupCount = getSampleGroups().size();
        final int margins = (groupCount - 1) * 3;
        squishedHeight = sampleCount == 0 || showGenotypes == false ? DEFAULT_SQUISHED_HEIGHT :
                Math.min(DEFAULT_SQUISHED_HEIGHT, Math.max(1, (getHeight() - getVariantBandHeight() - margins) / sampleCount));
        rowHeight = DEFAULT_EXPANDED_GENOTYPE_HEIGHT;

        // Set visibility window.  These values are appropriate for human dbsnp/1kg files, probably conservative otherwise
        // Ugly test on source is to avoid having to add "isIndexed" to a zillion feature source classes.  The intent
        // is to skip this if using a non-indexed source.
        if (!(source instanceof TribbleFeatureSource && ((TribbleFeatureSource) source).isIndexed() == false)) {
            int defVisibilityWindow = prefMgr.getAsInt(DEFAULT_VISIBILITY_WINDOW);
            if (defVisibilityWindow > 0) {
                setVisibilityWindow(defVisibilityWindow * 1000);
            } else {
                int vw = Math.max(10000, (100000 - 100 * (sampleCount() - 1)));
                setVisibilityWindow(vw);
            }
        }
    }

    @Override
    public int getHeight() {
        return !isVisible() ? 0 : height == 0 ? Math.min(DEFAULT_MAX_HEIGHT, getContentHeight()) : height;
    }

    @Override
    public TrackType getType() {
        return TrackType.variant;
    }

    @Override
    public boolean isFilterable() {
        return false;
    }

    private boolean defaultShowGenotypes() {
        return sampleCount() > 0;
    }

    public boolean isEnableMethylationRateSupport() {
        return enableMethylationRateSupport;
    }


    /**
     * Returns the height of a single sample (genotype) band
     *
     * @return
     */
    public int getGenotypeBandHeight() {
        switch (getDisplayMode()) {
            case SQUISHED:
                return getSquishedHeight();
            case COLLAPSED:
                return 0;
            default:
                return Math.max(1, Math.round(rowHeight));
        }
    }

    /**
     * Returns the total height of the track (including all sample/genotypes)
     *
     * @return
     */
    public int getContentHeight() {

        if (!isVisible()) {
            return 0;
        }

        int sampleCount = sampleCount();
        int h;
        if (getDisplayMode() == DisplayMode.COLLAPSED || sampleCount == 0 || showGenotypes == false) {
            h = getVariantsHeight();
        } else {
            final int groupCount = getSampleGroups().size();
            int margins = groupCount * 3;
            h = getVariantsHeight() + margins + Math.round(sampleCount * rowHeight);
        }
        return Math.max(WG_TRACK_HEIGHT, h);
    }

    public Object getHeader() {
        if (source instanceof TribbleFeatureSource) {
            return ((TribbleFeatureSource) source).getHeader();
        }
        return null;
    }

    /**
     * Return the height of the variant section only (no samples/genotypes)
     *
     * @return
     */
    private int getVariantsHeight() {
        return getVariantBandHeight() * getNumberOfFeatureLevels();
    }


    public void render(RenderContext context) {

        // Draw entire track.  TODO use clipBounds or visibleRect.
        Rectangle renderRect = context.getTrackRectangle();
        renderFeatures(context, renderRect);
    }

    protected void renderFeatures(RenderContext context, Rectangle ignore) {

        if (log.isTraceEnabled()) {
            String msg = String.format("renderFeatures: %s frame: %s", getName(), context.getReferenceFrame().getName());
            log.trace(msg);
        }

        PackedFeatures packedFeatures = packedFeaturesMap.get(context.getReferenceFrame());

        if (packedFeatures == null || !packedFeatures.overlapsInterval(context.getChr(), (int) context.getOrigin(), (int) context.getEndLocation() + 1)) {
            return;
        }

        try {
            renderFeatureImpl(context, packedFeatures);
        } catch (TribbleException e) {
            log.error("Tribble error", e);

            //Error loading features.  We'll let the user decide if this is "fatal" or not.
            boolean unload = MessageUtils.confirm("<html> Error loading features: " + e.getMessage() +
                    "<br>Unload track " + getName() + "?");
            if (unload) {
                Collection<Track> tmp = Arrays.asList((Track) this);
                IGV.getInstance().deleteTracks(tmp);
                IGV.getInstance().repaint();
            }

        }
    }


    /**
     * Render variants and genotypes if present.  Each site may contain multiple variants, which are rendered in
     * separate rows in expanded and squished modes.  Genotypes are rendered in a separate section below the variants
     * with each sample on a separaterow with alternating background.  The sample names are drawn in the name panel,
     * but the background color for each sample is drawn here.
     *
     * @param context
     * @param packedFeatures
     */
    @Override
    protected void renderFeatureImpl(RenderContext context, PackedFeatures packedFeatures) {

        Graphics2D g2D = context.getGraphics();

        Rectangle trackRectangle = context.trackRectangle;
        Rectangle clipBounds = new Rectangle(context.getClipBounds());

        Rectangle variantRect = new Rectangle(trackRectangle.x, trackRectangle.y, trackRectangle.width, getVariantsHeight());
        Rectangle genotypeRect = new Rectangle(trackRectangle.x, trackRectangle.y + getVariantsHeight(), trackRectangle.width, getGenotypeBandHeight());

        drawBackground(g2D, genotypeRect, clipBounds, BackgroundType.DATA);

        List<PackedFeatures.FeatureRow> rows = packedFeatures.getRows();

        if (rows.size() > 0) {

            final double locScale = context.getScale();
            final double origin = context.getOrigin();
            final double pXMin = variantRect.getMinX();
            final double pXMax = variantRect.getMaxX();

            for (PackedFeatures.FeatureRow row : rows) {

                List<VCFVariant> features = row.getFeatures();
                for (VCFVariant feature : features) {
                    VCFVariant variant = feature;

                    if (hideFiltered && variant.isFiltered()) {
                        continue;
                    }

                    int start = variant.getStart();
                    int end = variant.getEnd();
                    int pX = (int) ((start - origin) / locScale);
                    int dX = (int) Math.max(2, (end - start) / locScale);

                    if (pX + dX < pXMin) {
                        continue;
                    }
                    if (pX > pXMax) {
                        break;
                    }
                    int w = dX;
                    int x = pX;

                    if (w < 3) {
                        w = 3;
                        x--;
                    }

                    // if pixel width > 5 pixels create gap between variants
                    if (w > 5) {
                        x++;
                        w -= 2;
                    }

                    if (variantRect.intersects(clipBounds)) {
                        renderer.renderSiteBand(variant, variantRect, x, w, context);
                    }


                    if (showGenotypes) {

                        // Reset y position for each variant, it will be incremented as we loop through the samples
                        genotypeRect.y = trackRectangle.y + getVariantsHeight();
                        boolean hasGroups = getSampleGroups().size() > 1;


                        if (hasGroups && lastClipBounds != null) {
                            // Expand the clip bounds to be sure we clear previous labels, but not outside the bounds
                            // of the visible rectangle
                            Rectangle visibleRect = context.getVisibleRect();

                            int expandedTop = Math.max(visibleRect.y, Math.min(lastClipBounds.y, clipBounds.y));
                            int expandedBottom = Math.min(visibleRect.y + visibleRect.height,
                                    Math.max(clipBounds.y + clipBounds.height, lastClipBounds.y + lastClipBounds.height));
                            clipBounds.y = expandedTop;
                            clipBounds.height = expandedBottom - expandedTop;
                            context.getGraphics().setClip(clipBounds);
                        }

                        for (SampleGroup sampleGroup : getSampleGroups()) {
                            int yLabel = genotypeRect.y;
                            for (String sample : sampleGroup.samples()) {
                                if (genotypeRect.y > clipBounds.y + clipBounds.height) {
                                    break;
                                }
                                if (genotypeRect.y + genotypeRect.height > clipBounds.y) {
                                    renderer.renderGenotypeBandSNP(variant, context, genotypeRect, x, w, sample, genotypeColorMode, hideFiltered);
                                }
                                genotypeRect.y += genotypeRect.height;
                            }

                            String label = sampleGroup.label();
                            if (hasGroups && label != null) {
                                var r = new Rectangle(variantRect.x, yLabel, genotypeRect.width, Math.min(20, genotypeRect.y - yLabel));
                                GraphicUtils.drawVerticallyCenteredText(label, 0, r, context.getGraphics(), false, true);
                                drawGroupDivider(context.getGraphics(), genotypeRect, genotypeRect.y);
                            }

                            genotypeRect.y += groupGap;
                        }

                        boolean isSelected = selectedVariant != null && selectedVariant == variant;
                        if (isSelected) {
                            Graphics2D selectionGraphics = context.getGraphic2DForColor(Color.black);
                            selectionGraphics.drawRect(x, 0, w, this.getContentHeight());
                        }
                    }

                    lastClipBounds = context.getClipBounds();

                }
            }
        } else {
            g2D.setColor(Color.gray);
            GraphicUtils.drawCenteredText("No Variants Found", clipBounds, g2D);
        }

        renderBoundaryLines(g2D, clipBounds);

    }

    private void drawVariantBandBorder(Graphics2D g2D, Rectangle visibleRectangle, int variantBandY, int left, int right) {
        if (sampleCount() > 0 && showGenotypes && variantBandY >= visibleRectangle.y && variantBandY <= visibleRectangle.getMaxY()) {
            g2D.setColor(Color.lightGray);
            g2D.drawLine(left, variantBandY, right, variantBandY);
        }

    }


    /**
     * Render the border between variants / genotypes
     *
     * @param g2D
     * @param visibleRectangle
     */
    private void renderBoundaryLines(Graphics2D g2D, Rectangle visibleRectangle) {

        final int left = 0;
        final int right = visibleRectangle.width;

        // Variant / Genotype border
        if (sampleCount() > 0 && showGenotypes) {
            int variantGenotypeBorderY = getVariantsHeight();
            drawVariantBandBorder(g2D, visibleRectangle, variantGenotypeBorderY, left, right);
        }
    }

    /**
     * Render the name panel.
     * <p/>
     * NOTE:  The sample names are drawn in the drawBackground method.
     *
     * @param g
     * @param trackRectangle
     * @param visibleRect
     */
    @Override
    public void renderName(Graphics2D g, Rectangle trackRectangle, Rectangle visibleRect) {

        Graphics2D g2D = null;

        try {
            g2D = (Graphics2D) g.create();
            Rectangle clipRect = g2D.getClipBounds();

            g2D.setFont(FontManager.getFont(getFontSize()));
            g2D.setColor(Color.black);

            //   if(visibleRect.y < getVariantsHeight()) {
            Rectangle variantRect = new Rectangle(trackRectangle);
            variantRect.y = 0;
            variantRect.height = getVariantsHeight();
            GraphicUtils.drawWrappedText(getName(), variantRect, g2D, false);

            // The sample bounds list will get reset when  the names are drawn.
            sampleBounds.clear();
            drawBackground(g2D, trackRectangle, clipRect, BackgroundType.NAME);

            renderBoundaryLines(g2D, trackRectangle);
        } finally {
            if (g2D != null) {
                g2D.dispose();
            }
        }

    }

    @Override
    public int getSampleHeight() {
        return getGenotypeBandHeight();
    }

    @Override
    public int getNumRows() {
        return sampleCount();
    }

    @Override
    public int getSampleOffset() {
        return getVariantsHeight();
    }

    /**
     * Draws the "greenbar" type background.  Also draws the sample names.
     *
     * @param g2D
     * @param ignore
     * @param type
     */
    private void drawBackground(Graphics2D g2D, Rectangle trackRectangle, Rectangle ignore, BackgroundType type) {

        Rectangle clipBounds = g2D.getClipBounds();

        if (getDisplayMode() == DisplayMode.COLLAPSED || showGenotypes == false || trackRectangle.height < 2) {
            return;
        }

        // Create a rectangle for the genotype bands.  The height will be set to the band height, but the
        // y position will be incremented as we loop through the samples
        Rectangle bandRectangle = new Rectangle(trackRectangle);
        bandRectangle.y = getVariantsHeight(); // Start below the variant bands
        bandRectangle.height = getGenotypeBandHeight();

        int bandFontSize = Math.min(getFontSize(), (int) bandRectangle.getHeight() - 1);
        Font font = FontManager.getFont(bandFontSize);
        Font oldFont = g2D.getFont();
        g2D.setFont(font);

        boolean supressFill = trackRectangle.height < 4;
        boolean hasGroups = getSampleGroups().size() > 1;
        boolean b = false;

        for (SampleGroup sampleGroup : getSampleGroups()) {

            for (String sample : sampleGroup.samples()) {

                if (bandRectangle.y > clipBounds.y + clipBounds.height) {
                    break;
                }

                Color bgColor = b ? BAND1_COLOR : BAND2_COLOR;
                b = !b;

                if (bandRectangle.y + bandRectangle.height > clipBounds.y) {

                    g2D.setColor(bgColor); //darkMode ? UIManager.getColor("Panel.background") : Color.white);

                    if (!supressFill) {
                        g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, bandRectangle.height);
                    }

                    if (type == BackgroundType.NAME) {

                        sampleBounds.add(new SampleBounds(bandRectangle.y, bandRectangle.y + bandRectangle.height, sample));

                        if (bandRectangle.height >= 3) {
                            String printName = sample;
                            g2D.setColor(darkMode ? Color.white : Color.black);
                            GraphicUtils.drawWrappedText(printName, bandRectangle, g2D, false);
                        }
                    }
                }

                bandRectangle.y += bandRectangle.height;
            }
            if (hasGroups) {
                drawGroupDivider(g2D, trackRectangle, bandRectangle.y);
            }
            bandRectangle.y += groupGap;
        }
        g2D.setFont(oldFont);
    }


    public boolean getHideFiltered() {
        return hideFiltered;
    }

    public void setHideFiltered(boolean value) {
        this.hideFiltered = value;
    }

    public ColorMode getGenotypeColorMode() {
        return genotypeColorMode;
    }

    public void setGenotypeColorMode(ColorMode mode) {
        this.genotypeColorMode = mode;
    }

    public ColorMode getSiteColorMode() {
        return siteColorMode;
    }

    public void setSiteColorMode(ColorMode siteColorMode) {
        this.siteColorMode = siteColorMode;
    }

    @Override
    public void setColor(Color color) {
        // Setting color implicitly turns of "color by" modes
        this.genotypeColorMode = ColorMode.NONE;
        this.siteColorMode = ColorMode.NONE;
        super.setColor(color);
    }

    public String getTooltipText(int y) {
        if (y < getVariantsHeight()) {
            return super.getTooltipText(y);
        } else {
            String sample = getSampleAtPosition(y);
            return sample;
        }
    }

    /**
     * Handle a mouse click from the data panel.
     *
     * @param te - wraps the MouseClickEvent and reference frame.
     * @return true if the click is handled, false otherwise
     */
    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        if (te.getMouseEvent().isPopupTrigger()) {
            return false;
        }
        if (IGV.getInstance().isShowDetailsOnClick()) {
            openTooltipWindow(te);
        }

        return true;
    }

    /**
     * Return popup text for the given position
     *
     * @param chr
     * @param position - position in UCSC "0 based"  genomic coordinates
     * @param mouseX
     * @param frame    @return
     */
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        try {
            double maxDistance = 10 * frame.getScale();
            if (mouseY < getVariantsHeight()) {
                int modY = mouseY;
                Variant variant = getFeatureClosest(position, modY, frame, maxDistance);
                if (variant == null) return null;
                return getVariantToolTip(variant);
            } else {
                if (sampleBounds == null || sampleBounds.isEmpty()) return null;
                String sample = getSampleAtPosition(mouseY);
                if (sample == null) return null;

                Variant variant = getFeatureClosest(position, -1, frame, maxDistance);
                return getSampleToolTip(sample, variant);
            }
        } catch (Exception e) {
            log.error("Error getting value string", e);
            return null;
        }
    }

    /**
     * Return the sample at the give pixel position
     *
     * @param y - screen position in pixels
     * @return
     */
    private String getSampleAtPosition(int y) {

        if (sampleBounds.isEmpty()) {
            return null;
        }
        String sample = null;

        // Estimate the index of the sample, then do a linear search
        final int sampleCount = sampleBounds.size();

        int firstSampleY = sampleBounds.get(0).top;
        int idx = Math.max(0, Math.min((y - firstSampleY) / getGenotypeBandHeight(), sampleCount - 1));

        SampleBounds bounds = sampleBounds.get(idx);
        if (bounds.contains(y)) {
            sample = bounds.sample;
        } else if (bounds.top > y) {
            while (idx > 0) {
                idx--;
                bounds = sampleBounds.get(idx);
                if (bounds.contains(y)) {
                    sample = bounds.sample;
                }
            }
        } else {
            while (idx < sampleCount - 1) {
                idx++;
                bounds = sampleBounds.get(idx);
                if (bounds.contains(y)) {
                    sample = bounds.sample;
                }
            }
        }
        return sample;
    }

    /**
     * Return the variant closest to the genomic position in the given reference frame, within the prescribed tolerance
     *
     * @param position
     * @param y           pixel position in panel coordinates (i.e. not track coordinates)
     * @param frame
     * @param maxDistance
     * @return
     */
    protected Variant getFeatureClosest(double position, int y, ReferenceFrame frame, double maxDistance) {

        PackedFeatures<PackedFeature> packedFeatures = packedFeaturesMap.get(frame);

        if (packedFeatures == null) {
            return null;
        }

        Feature feature = null;
        List<? extends Feature> features;

        //We search only the specified row if y is a meaningful value.
        //Otherwise we search everything
        int row = (y / getVariantBandHeight());
        if (y < 0 || row >= getNumberOfFeatureLevels()) {
            features = packedFeatures.getFeatures();
        } else {
            features = packedFeatures.getRows().get(row).getFeatures();
        }

        if (features != null) {
            feature = FeatureUtils.getFeatureClosest(position, features);
        }
        if (feature == null ||
                ((position < feature.getStart() - maxDistance) || (position > feature.getEnd() + maxDistance))) {
            return null;
        } else {
            return (Variant) feature;
        }


    }

    private String getVariantToolTip(Variant variant) {
        String id = variant.getID();

        StringBuffer toolTip = new StringBuffer();
        if (id.length() > 0) {
            toolTip.append("ID: " + id + "<br>");
        }
        toolTip.append("Chr: " + variant.getChr());
        toolTip.append("<br>Position: " + variant.getPositionString());
        toolTip.append("<br>Reference: " + variant.getReference());
        List<Allele> alternates = variant.getAlternateAlleles();
        String alternateString = null;
        if (alternates.size() > 0) {
            String tmp = alternates.get(0).toString();
            alternateString = StringUtils.join(alternates, ",");
            toolTip.append("<br>Alternate: " + alternateString);
        }

        double qual = variant.getPhredScaledQual();
        String qualString = variant.hasLog10PError() ? numFormat.format(qual) : ".";
        toolTip.append("<br>Qual: " + qualString);
        toolTip.append("<br>Type: " + variant.getType());
        if (variant.isFiltered()) {
            toolTip.append("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.append(getFilterTooltip(variant));
        } else {
            toolTip.append("<br>Is Filtered Out: No</b><br>");
        }

        if (alternateString != null) {
            toolTip.append("<br><b>Alleles:</b>");

            toolTip.append("<br>Alternate Alleles: " + alternateString);

            int[] ac = variant.getAlleleCounts();
            if (ac != null) {
                String acString = ac.length > 1 ? "<br>Allele Counts: " : "<br>Allele Count: ";
                for (int i = 0; i < ac.length; i++) {
                    acString += Integer.toString(ac[i]);
                    if (i < ac.length - 1) acString += ", ";
                }
                toolTip.append(acString);
            }

            int totalAlleleCount = variant.getTotalAlleleCount();
            if (totalAlleleCount > 0) {
                toolTip.append("<br>Total # Alleles: " + String.valueOf(totalAlleleCount));
            }

            double[] af = variant.getAlleleFreqs();

            int nonNegativeCounts = 0;
            for (int i = 0; i < af.length; i++) {
                if (af[i] >= 0) nonNegativeCounts++;
            }
            if (nonNegativeCounts > 0) {
                String afString = nonNegativeCounts > 1 ? "<br>Allele Fequencies: " : "<br>Allele Frequency: ";
                for (int i = 0; i < af.length; i++) {
                    if (af[i] >= 0) {
                        afString += Double.toString(af[i]);
                        if (i < af.length - 1) afString += ", ";
                    }
                }
                toolTip.append(afString);
            }
        }
        if (variant.getAttributes().size() > 0) {
            toolTip.append(getVariantInfo(variant));
        }


        return toolTip.toString();
    }

    protected String getVariantInfo(Variant variant) {
        Set<String> keys = variant.getAttributes().keySet();
        if (keys.size() > 0) {
            String toolTip = "<br><br><b>Variant Attributes</b>";
            int count = 0;

            // Put AF and GMAF and put at the 0, if present
            String k = "AF";
            String afValue = variant.getAttributeAsString(k);
            if (afValue != null && afValue.length() > 0 && !afValue.equals("null")) {
                toolTip = toolTip.concat("<br>" + getFullName(k) + ": " + variant.getAttributeAsString(k));
            }

            k = "GMAF";
            afValue = variant.getAttributeAsString(k);
            if (afValue != null && afValue.length() > 0 && !afValue.equals("null")) {
                toolTip = toolTip.concat("<br>" + getFullName(k) + ": " + variant.getAttributeAsString(k));
            }
            int maxFilterLines = getMaxFilterLines();
            for (String key : keys) {
                count++;

                if (key.equals("AF") || key.equals("GMAF")) continue;

                if (count > maxFilterLines) {
                    toolTip = toolTip.concat("<br>....");
                    break;
                }
                toolTip = toolTip.concat("<br>" + getFullName(key) + ": " + variant.getAttributeAsString(key));

            }
            return toolTip;
        }
        return " ";
    }

    /**
     * The maximum number of filter lines to show for variants.
     * We show more info if the user is displaying a separate window than
     * if using tooltip
     *
     * @return
     */
    private int getMaxFilterLines() {
        return IGV.getInstance().isShowDetailsOnHover() ? MAX_FILTER_LINES : 1000;
    }

    private String getGenotypeInfo(Genotype genotype) {
        final Map<String, Object> attributes = genotype.getAttributes();
        Set<String> keys = attributes.keySet();
        if (keys.size() > 0) {
            String tooltip = "<br><b>Genotype Attributes</b>";
            for (String key : keys) {
                tooltip = tooltip.concat("<br>" + getFullName(key) + ": " + attributes.get(key));
            }
            return tooltip;
        }
        return null;
    }

    public void clearSelectedVariant() {
        selectedVariant = null;
    }

    public int getSquishedHeight() {
        return squishedHeight;
    }

    public void setSquishedHeight(int squishedHeight) {
        this.squishedHeight = squishedHeight;
    }


    public boolean isShowGenotypes() {
        return showGenotypes;
    }

    public void setShowGenotypes(boolean showGenotypes) {
        this.showGenotypes = showGenotypes;
    }

    /**
     * The height of the top band representing the variant call
     */
    public int getVariantBandHeight() {
        return getDisplayMode() == DisplayMode.SQUISHED ? DEFAULT_SQUISHED_VARIANT_HEIGHT : DEFAULT_EXPANDED_VARIANT_HEIGHT;
    }

    public enum ColorMode {
        GENOTYPE, METHYLATION_RATE, ALLELE_FREQUENCY, NONE, ALLELE_FRACTION
    }

    public static enum BackgroundType {
        NAME, ATTRIBUTE, DATA;
    }


    static Map<String, String> fullNames = new HashMap();

    static {
        fullNames.put("AA", "Ancestral Allele");
        fullNames.put("AC", "Allele Count");
        fullNames.put("AN", "Total Alleles");
        fullNames.put("AF", "Allele Frequency");
        fullNames.put("DP", "Depth");
        fullNames.put("MQ", "Mapping Quality");
        fullNames.put("NS", "Number of Samples with Data");
        fullNames.put("BQ", "RMS Base Quality");
        fullNames.put("SB", "Strand Bias");
        fullNames.put("DB", "dbSNP Membership");
        fullNames.put("GQ", "Genotype Quality");
        fullNames.put("GL", "Genotype Likelihoods");  //Hom-ref, het, hom-var break this down into a group
    }

    static String getFullName(String key) {
        return fullNames.containsKey(key) ? fullNames.get(key) : key;
    }


    private String getSampleToolTip(String sample, Variant variant) {

        if (variant == null) return null;
        double goodBaseCount = variant.getGenotype(sample).getAttributeAsDouble("GB");
        if (Double.isNaN(goodBaseCount)) goodBaseCount = 0;
        if (isEnableMethylationRateSupport() && goodBaseCount < 10) {
            return sample;
        }
        String id = variant.getID();
        StringBuffer toolTip = new StringBuffer();
        toolTip = toolTip.append("Chr: " + variant.getChr());
        toolTip = toolTip.append("<br>Position: " + variant.getPositionString());
        toolTip = toolTip.append("<br>ID: " + id + "<br>");
        toolTip = toolTip.append("<br><b>Genotype Information</b>");
        toolTip = toolTip.append("<br>Sample: " + sample);

        Genotype genotype = variant.getGenotype(sample);
        if (genotype != null) {
            toolTip = toolTip.append("<br>Genotype: " + genotype.getGenotypeString());
            toolTip = toolTip.append("<br>Quality: " + numFormat.format(genotype.getPhredScaledQual()));
            toolTip = toolTip.append("<br>Type: " + genotype.getTypeString());
        }
        if (variant.isFiltered()) {
            toolTip = toolTip.append("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.append(getFilterTooltip(variant));
        } else {
            toolTip = toolTip.append("<br>Is Filtered Out: No</b><br>");
        }

        if (genotype != null) {
            String sInfoStr = getGenotypeInfo(genotype);
            if (sInfoStr != null) {
                toolTip = toolTip.append(sInfoStr + "<br>");
            }
        }
        return toolTip.toString();
    }


    private String getFilterTooltip(Variant variant) {
        Collection filters = variant.getFilters();
        String toolTip = "<br>";
        for (Object filter : filters) {
            toolTip = toolTip.concat("- " + (String) filter + "<br>");
        }

        return toolTip;
    }


    /**
     * Return the {@code Variant} object closest to the specified event
     *
     * @param te
     * @return
     * @api
     */
    public Variant getSelectedVariant(final TrackClickEvent te) {
        final ReferenceFrame referenceFrame = te.getFrame();
        Variant selVariant = null;
        if (referenceFrame != null && referenceFrame.getName() != null) {
            final double position = te.getChromosomePosition();
            double maxDistance = 10 * referenceFrame.getScale();
            selVariant = getFeatureClosest(position, te.getMouseEvent().getY(), referenceFrame, maxDistance);
        }
        return selVariant;
    }


    @Override
    public List<Component> getPopupMenuItems(final TrackClickEvent te) {
        selectedVariant = getSelectedVariant(te);
        if (selectedVariant != null) {
            repaint();
        }
        return VariantTrackMenuHelper.getMenuItems(this, selectedVariant, te);
    }

    @Override
    public boolean hasDisplayMode() {
        return true;
    }

    /**
     * Return the index for the sample.  This is a very inefficient implementation, but we don't care because
     * these lists are tiny.
     *
     * @param sample
     * @return
     */
    private int getSampleIndex(String sample) {
        for (int i = 0; i < sampleBounds.size(); i++) {
            if (sample.equals(sampleBounds.get(i).sample)) {
                return i;
            }
        }
        return -1;
    }


    /**
     * Return the nextLine or previous feature relative to the center location.
     * <p/>
     * Loop through "next feature from super implementation until first non-filtered variant is found.
     *
     * @param chr
     * @param center
     * @param forward
     * @return
     * @throws IOException
     */
    @Override
    public Feature nextFeature(String chr, double center, boolean forward, ReferenceFrame frame) throws IOException {

        if (getHideFiltered()) {
            Feature f;
            while ((f = super.nextFeature(chr, center, forward, frame)) != null) {
                if (!(f instanceof Variant) || !((Variant) f).isFiltered()) {
                    return f;
                } else {
                    chr = f.getChr();
                    center = (f.getStart() + f.getEnd()) / 2 + 1;
                }
            }
            return null;
        } else {
            return super.nextFeature(chr, center, forward, frame);
        }
    }

    static class SampleBounds {
        int top;
        int bottom;
        String sample;

        SampleBounds(int top, int bottom, String sample) {
            this.top = top;
            this.bottom = bottom;
            this.sample = sample;
        }

        boolean contains(int y) {
            return y >= top && y <= bottom;
        }
    }

    void sendToCircularView(TrackClickEvent e) {

        List<Feature> visibleFeatures;
        if (e.getFrame() == null) {
            visibleFeatures = new ArrayList<>();
            for (ReferenceFrame frame : FrameManager.getFrames()) {
                visibleFeatures.addAll(getVisibleFeatures(frame));
            }
        } else {
            visibleFeatures = getVisibleFeatures(e.getFrame());
        }

        List<Feature> svFeatures = visibleFeatures.stream().filter(f -> {
            Variant v = f instanceof MateVariant ? ((MateVariant) f).mate : (Variant) f;
            Map<String, Object> attrs = v.getAttributes();
            return attrs.containsKey("CHR2") && attrs.containsKey("END");
        }).collect(Collectors.toList());

        if (svFeatures.isEmpty()) {
            MessageUtils.showMessage("No structural variants found.");
        } else {
            CircularViewUtilities.sendVariantsToJBrowse(svFeatures, getName(), CIRC_VIEW_DEFAULT_COLOR);
        }
    }

    @Override
    public List<Feature> getVisibleFeatures(ReferenceFrame frame) {
        if (frame.getChrName().equals(Globals.CHR_ALL) &&
                this.source instanceof TribbleFeatureSource.NonIndexedFeatureSource) {
            try {
                return ((TribbleFeatureSource.NonIndexedFeatureSource) this.source).getAllFeatures();
            } catch (IOException e) {
                return Collections.emptyList();
            }
        } else {
            return super.getVisibleFeatures(frame);
        }
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("showGenotypes")) {
            this.showGenotypes = Boolean.parseBoolean(element.getAttribute("showGenotypes"));
        }

        if (element.hasAttribute("squishedHeight")) {
            this.squishedHeight = Integer.parseInt(element.getAttribute("squishedHeight"));
        }

        if (element.hasAttribute("genotypeColorMode")) {
            this.genotypeColorMode = ColorMode.valueOf(element.getAttribute("genotypeColorMode"));
        } else if (element.hasAttribute("coloring")) {
            // backward compatibility
            this.genotypeColorMode = ColorMode.valueOf(element.getAttribute("coloring"));
        }

        if (element.hasAttribute("siteColorMode")) {
            this.siteColorMode = ColorMode.valueOf(element.getAttribute("siteColorMode"));

        }
    }


    @Override
    public void marshalJSON(org.json.JSONObject json) {

        super.marshalJSON(json);

        if (showGenotypes != defaultShowGenotypes()) {
            json.put("showGenotypes", showGenotypes);
        }
        if (this.squishedHeight != DEFAULT_SQUISHED_HEIGHT) {
            json.put("squishedHeight", squishedHeight);
        }
        if (genotypeColorMode != ColorMode.GENOTYPE) {
            json.put("genotypeColorMode", genotypeColorMode.toString());
        }
        if (siteColorMode != null) {
            json.put("siteColorMode", siteColorMode.toString());
        }
    }


    @Override
    public void unmarshalJSON(org.json.JSONObject json) {

        super.unmarshalJSON(json);

        if (json.has("showGenotypes")) {
            this.showGenotypes = json.getBoolean("showGenotypes");
        }

        if (json.has("squishedHeight")) {
            this.squishedHeight = json.getInt("squishedHeight");
        }

        if (json.has("genotypeColorMode")) {
            this.genotypeColorMode = ColorMode.valueOf(json.getString("genotypeColorMode"));
        }

        if (json.has("siteColorMode")) {
            this.siteColorMode = ColorMode.valueOf(json.getString("siteColorMode"));
        }
    }

}
