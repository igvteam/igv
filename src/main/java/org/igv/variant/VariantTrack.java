//chr2:128,565,093-128,565,156

package org.igv.variant;

import htsjdk.tribble.Feature;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.igv.Globals;
import org.igv.event.IGVEventObserver;
import org.igv.feature.FeatureUtils;
import org.igv.feature.PackedFeature;
import org.igv.circview.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.GraphicUtils;
import org.igv.sample.SampleGroup;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.color.ColorPalette;
import org.igv.ui.color.ColorUtilities;
import org.igv.ui.color.PaletteColorTable;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.StringUtils;
import org.igv.variant.vcf.MateVariant;
import org.igv.variant.vcf.VCFVariant;
import org.w3c.dom.Element;

import javax.swing.UIManager;
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
    // Alternating "greenbar" background for the genotype bands.  In dark mode these are derived from the panel
    // background so the bands stay dark -- the sample names drawn over them are white.
    private static final Color BAND1_COLOR;
    private static final Color BAND2_COLOR;

    static {
        if (Globals.isDarkMode()) {
            Color background = UIManager.getColor("Panel.background");
            if (background == null) background = new Color(60, 63, 65);
            BAND1_COLOR = ColorUtilities.shiftBrightness(background, 12);
            BAND2_COLOR = background;
        } else {
            BAND1_COLOR = new Color(245, 245, 245);
            BAND2_COLOR = Color.white;
        }
    }

    private final static int DEFAULT_EXPANDED_GENOTYPE_HEIGHT = 15;
    private final static int VARIANT_BAND_HEIGHT = 25;
    private final static int MAX_FILTER_LINES = 15;
    private final static int WG_TRACK_HEIGHT = 40;
    private final static int DEFAULT_SQUISHED_GENOTYPE_HEIGHT = 4;

    /**
     * Color for variants with no value for the INFO attribute being colored by.
     */
    private static final Color NO_ATTRIBUTE_VALUE_COLOR = Color.gray;

    /**
     * Palette for attribute values no color scheme covers.  Matches the default used by igv.js.
     */
    private static final String ATTRIBUTE_PALETTE = "Set 1";

    /**
     * How far apart, in RGB, two attribute colors have to be to count as distinguishable.
     */
    private static final int MIN_COLOR_DISTANCE = 60;

    /**
     * How many generated colors to try before settling for the furthest one found.
     */
    private static final int MAX_COLOR_ATTEMPTS = 500;

    /**
     * INFO attribute types that can be colored by.
     */
    private static final Set<VCFHeaderLineType> COLORABLE_INFO_TYPES = EnumSet.of(
            VCFHeaderLineType.String,
            VCFHeaderLineType.Character,
            VCFHeaderLineType.Integer,
            VCFHeaderLineType.Float,
            VCFHeaderLineType.Flag);

    /**
     * INFO attribute types that are quantities rather than categories.  Coloring these by value would give a
     * color per variant, so a color scale is defined the first time one is selected.
     */
    private static final Set<VCFHeaderLineType> NUMERIC_INFO_TYPES = EnumSet.of(
            VCFHeaderLineType.Integer,
            VCFHeaderLineType.Float);


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
     * Current coloring option
     */
    private ColorMode genotypeColorMode = ColorMode.GENOTYPE;


    private ColorMode siteColorMode;

    /**
     * The VCF INFO attribute the variant band is colored by.  Only used when siteColorMode == ATTRIBUTE.
     */
    private String colorByAttribute;

    /**
     * Color tables for "color by INFO attribute", keyed by attribute ID.  Tables are created on demand and
     * accumulate assignments as attribute values are encountered, so they are per-track state.
     */
    private final Map<String, PaletteColorTable> attributeColorTables = new HashMap<>();

    /**
     * Colors the user chose explicitly for attribute values, keyed by attribute ID then by value.  These take
     * precedence over color schemes, which in turn take precedence over colors assigned from the palette.
     */
    private final Map<String, Map<String, Color>> attributeColorOverrides = new HashMap<>();

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

        this.defaultExpandedRowHeight = DEFAULT_EXPANDED_GENOTYPE_HEIGHT;
        this.defaultSquishedRowHeight = DEFAULT_SQUISHED_GENOTYPE_HEIGHT;
        setDisplayMode(DisplayMode.EXPANDED);

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
    public void minimizeHeight() {
        setRowHeight(1);
        int newHeight = Math.max(getContentHeight(), getMinimumHeight());
        setHeight(Math.min(newHeight, getHeight()));
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
        return Math.max(1, rowHeight);
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
            h = getVariantsHeight() + margins + sampleCount * getGenotypeBandHeight();
        }
        return Math.max(WG_TRACK_HEIGHT, h);
    }

    public Object getHeader() {
        try {
            return source == null ? null : source.getHeader();
        } catch (IOException e) {
            log.error("Error fetching header for " + getName(), e);
            return null;
        }
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

                List<Variant> features = row.getFeatures();
                for (Variant feature : features) {
                    Variant variant = feature;

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
                                context.getGraphics().setColor(UIConstants.getTrackPanelForeground());
                                GraphicUtils.drawVerticallyCenteredText(label, 0, r, context.getGraphics(), false, true);
                                drawGroupDivider(context.getGraphics(), genotypeRect, genotypeRect.y);
                            }

                            genotypeRect.y += groupGap;
                        }

                        boolean isSelected = selectedVariant != null && selectedVariant == variant;
                        if (isSelected) {
                            Graphics2D selectionGraphics =
                                    context.getGraphic2DForColor(UIConstants.getTrackPanelForeground());
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
            g2D.setColor(UIConstants.getTrackPanelForeground());

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
    public boolean hasRows() {
        return true;
    }

    @Override
    public int getReservedHeight() {
        return getVariantsHeight();
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

    /**
     * @return the VCF INFO attribute the variant band is colored by, or null if not coloring by attribute.
     */
    public String getColorByAttribute() {
        return colorByAttribute;
    }

    /**
     * Color the variant band by the values of a VCF INFO attribute.  Passing null reverts to a flat track color.
     *
     * @param key an INFO attribute ID, or null
     */
    public void setColorByAttribute(String key) {
        this.colorByAttribute = key;
        if (key == null) {
            if (siteColorMode == ColorMode.ATTRIBUTE) {
                this.siteColorMode = ColorMode.NONE;
            }
        } else {
            this.siteColorMode = ColorMode.ATTRIBUTE;
            assignColorsForLoadedFeatures(key);
        }
    }

    /**
     * Return the color for a variant when coloring by an INFO attribute.  Variants with no value for the
     * attribute are drawn gray, following igv.js.
     */
    public Color getAttributeColor(Variant variant) {

        if (colorByAttribute == null) {
            return getColor();
        }

        String value = normalizeAttributeValue(variant.getAttributeAsString(colorByAttribute));
        if (value == null) {
            return NO_ATTRIBUTE_VALUE_COLOR;
        }
        return getAttributeColor(colorByAttribute, value);
    }

    /**
     * Return the color for a value of an INFO attribute.  Colors the user chose for this track win over colors
     * from a scheme, which win over colors assigned from the palette -- so importing a scheme takes effect on
     * tracks that are already open, without discarding the user's own choices.
     */
    public Color getAttributeColor(String key, String value) {

        Color override = getAttributeColorOverride(key, value);
        if (override != null) {
            return override;
        }

        AbstractColorScale scale = VariantColorSchemes.getScale(key);
        if (scale != null) {
            Double number = numericValue(key, value);
            return number == null ? NO_ATTRIBUTE_VALUE_COLOR : scale.getColor(number.floatValue());
        }

        Color color = VariantColorSchemes.getColor(key, value);
        return color != null ? color : assignPaletteColor(key, value);
    }

    /**
     * Return the palette color for a value no scheme covers, assigning one if this is the first time it is seen.
     * <p>
     * Palette entries close to a color a scheme already uses for this attribute are skipped -- otherwise a value
     * with no scheme entry can come out looking like one that has one, which is worse than an unrelated color.
     */
    private Color assignPaletteColor(String key, String value) {

        PaletteColorTable colorTable = getAttributeColorTable(key);

        if (!colorTable.getColorMap().containsKey(value.toLowerCase())) {
            // Always store a color chosen here.  Leaving it to the table would undo the filtering below, as it
            // takes palette[size] regardless of whether that entry was skipped as too similar.
            colorTable.put(value, nextDistinctColor(key, colorTable));
        }
        return colorTable.get(value);
    }

    /**
     * A color that is neither already assigned for this attribute nor close to a color a scheme uses for it.
     * Palette colors are preferred; once they are used up or rejected, colors are generated until one is far
     * enough from everything in use.  Generation is deterministic, so the assignments are reproducible.
     */
    private Color nextDistinctColor(String key, PaletteColorTable colorTable) {

        List<Color> used = new ArrayList<>(colorTable.getColorMap().values());
        used.addAll(VariantColorSchemes.getColors(key));

        ColorPalette palette = ColorUtilities.getPalette(ATTRIBUTE_PALETTE);
        if (palette != null) {
            for (Color candidate : palette.getColors()) {
                if (minDistance(candidate, used) >= MIN_COLOR_DISTANCE) {
                    return candidate;
                }
            }
        }

        // Keep the furthest candidate seen, so that even a crowded attribute never repeats a color outright
        Color best = ColorUtilities.randomColor(used.size());
        double bestDistance = minDistance(best, used);

        for (int i = 1; i <= MAX_COLOR_ATTEMPTS && bestDistance < MIN_COLOR_DISTANCE; i++) {
            Color candidate = ColorUtilities.randomColor(used.size() + i);
            double distance = minDistance(candidate, used);
            if (distance > bestDistance) {
                best = candidate;
                bestDistance = distance;
            }
        }
        return best;
    }

    /**
     * Distance from a color to the nearest of those already in use, or a large value if none are.
     */
    private static double minDistance(Color color, List<Color> used) {
        double min = Double.MAX_VALUE;
        for (Color c : used) {
            min = Math.min(min, distance(color, c));
        }
        return min;
    }

    private static double distance(Color c1, Color c2) {
        int dr = c1.getRed() - c2.getRed();
        int dg = c1.getGreen() - c2.getGreen();
        int db = c1.getBlue() - c2.getBlue();
        return Math.sqrt(dr * dr + dg * dg + db * db);
    }

    /**
     * Set the color for one value of an INFO attribute on this track.  A null color clears the choice, reverting
     * to the scheme or palette color.
     */
    public void setAttributeColorOverride(String key, String value, Color color) {
        synchronized (attributeColorOverrides) {
            Map<String, Color> overrides = attributeColorOverrides.computeIfAbsent(key, k -> new LinkedHashMap<>());
            if (color == null) {
                overrides.remove(value.toLowerCase());
            } else {
                overrides.put(value.toLowerCase(), color);
            }
        }
    }

    /**
     * Discard the colors chosen for an INFO attribute on this track, reverting to scheme and palette colors.
     */
    public void clearAttributeColorOverrides(String key) {
        synchronized (attributeColorOverrides) {
            attributeColorOverrides.remove(key);
        }
    }

    /**
     * @return the colors chosen for values of an INFO attribute, keyed by lower case value.  Never null.
     */
    public Map<String, Color> getAttributeColorOverrides(String key) {
        synchronized (attributeColorOverrides) {
            Map<String, Color> overrides = attributeColorOverrides.get(key);
            return overrides == null ? Collections.emptyMap() : new LinkedHashMap<>(overrides);
        }
    }

    /**
     * Look up one chosen color.  This is on the per variant render path, so it reads the map under its lock
     * rather than copying it -- copying makes a repaint cost variants x overrides once anything is edited.
     */
    private Color getAttributeColorOverride(String key, String value) {
        synchronized (attributeColorOverrides) {
            Map<String, Color> overrides = attributeColorOverrides.get(key);
            return overrides == null ? null : overrides.get(value.toLowerCase());
        }
    }

    /**
     * @return true if the attribute is a quantity rather than a category, and so needs a color scale.
     */
    public boolean isNumericAttribute(String key) {
        Object header = getHeader();
        if (!(header instanceof VCFHeader)) {
            return false;
        }
        VCFInfoHeaderLine line = ((VCFHeader) header).getInfoHeaderLine(key);
        return line != null && NUMERIC_INFO_TYPES.contains(line.getType());
    }

    /**
     * Return the range of an attribute's values among the currently loaded features, as {min, max}, or null if
     * there are no numeric values.  Used to prefill the color scale editor with a range that suits the data.
     */
    public double[] getAttributeRange(String key) {

        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        for (String value : getAttributeValues(key)) {
            Double number = numericValue(key, value);
            if (number != null) {
                min = Math.min(min, number);
                max = Math.max(max, number);
            }
        }

        return min > max ? null : new double[]{min, max};
    }

    /**
     * Return the distinct values of an INFO attribute among the currently loaded features, sorted.  This is what
     * the color legend shows -- the values the user can actually see.
     */
    public SortedSet<String> getAttributeValues(String key) {

        SortedSet<String> values = new TreeSet<>(String.CASE_INSENSITIVE_ORDER);
        synchronized (packedFeaturesMap) {
            for (PackedFeatures<PackedFeature> packedFeatures : packedFeaturesMap.values()) {
                for (PackedFeature feature : packedFeatures.getFeatures()) {
                    if (feature instanceof Variant) {
                        String value = normalizeAttributeValue(((Variant) feature).getAttributeAsString(key));
                        if (value != null) {
                            values.add(value);
                        }
                    }
                }
            }
        }
        return values;
    }

    /**
     * Return the table of palette colors assigned to values of the given INFO attribute.  It holds only values no
     * color scheme covers -- those are assigned a color as they are encountered.
     */
    public PaletteColorTable getAttributeColorTable(String key) {
        synchronized (attributeColorTables) {
            return attributeColorTables.computeIfAbsent(key,
                    k -> new PaletteColorTable(ColorUtilities.getPalette(ATTRIBUTE_PALETTE)));
        }
    }

    /**
     * Assign palette colors to the attribute values of the currently loaded features, in sorted order.  Colors
     * are otherwise assigned in the order values happen to be drawn, which makes them depend on where the user
     * navigated first.  This does not affect values with a predefined color, and values encountered later are
     * appended as usual.
     */
    private void assignColorsForLoadedFeatures(String key) {

        if (VariantColorSchemes.getScale(key) != null) {
            return;     // Colors come from a continuous scale, nothing to assign
        }

        SortedSet<String> values = new TreeSet<>(String.CASE_INSENSITIVE_ORDER);

        synchronized (packedFeaturesMap) {
            for (PackedFeatures<PackedFeature> packedFeatures : packedFeaturesMap.values()) {
                for (PackedFeature feature : packedFeatures.getFeatures()) {
                    if (feature instanceof Variant) {
                        String value = normalizeAttributeValue(((Variant) feature).getAttributeAsString(key));
                        if (value != null && VariantColorSchemes.getColor(key, value) == null) {
                            values.add(value);
                        }
                    }
                }
            }
        }

        values.forEach(value -> assignPaletteColor(key, value));
    }

    /**
     * The number to color a variant by, for an attribute colored by a scale.  Null if the value holds no number.
     * <p>
     * A "Number=A" or "Number=R" attribute carries one value per allele, so a multi-allelic record arrives here
     * as a comma separated list.  Which values count, and how they combine, follows the header: the reference
     * allele's value in a Number=R list is skipped, and the alternate values are summed for an allele frequency
     * and otherwise reduced to their maximum (see {@link Aggregation}).
     */
    private Double numericValue(String key, String value) {

        if (value == null) {
            return null;
        }

        String[] parts = value.split(",");

        // Number=R lists the reference allele first.  Per-allele coloring is about the alternates -- a sum that
        // included the reference would not be the frequency of the non-reference alleles it claims to be.
        int first = parts.length > 1 && getCountType(key) == VCFHeaderLineCount.R ? 1 : 0;

        return aggregate(parts, first, getAggregation(key));
    }

    /**
     * How the per-allele values of an attribute combine into the one number a variant is colored by.
     */
    enum Aggregation {
        /**
         * The largest value.  The default: for a score it is the most severe allele, and it never runs off the
         * end of the scale.
         */
        MAX,
        /**
         * The total.  For an allele frequency the total across alternate alleles is what the attribute means,
         * and it matches the "Allele Frequency" color mode (see
         * {@link org.igv.variant.vcf.VCFVariant#getAlternateAlleleFrequency}).
         */
        SUM
    }

    private static Aggregation getAggregation(String key) {
        for (String frequencyKey : VCFVariant.ALLELE_FREQUENCY_KEYS) {
            if (frequencyKey.equalsIgnoreCase(key)) {
                return Aggregation.SUM;
            }
        }
        return Aggregation.MAX;
    }

    /**
     * The header's declared cardinality for an attribute, or UNBOUNDED if it is not declared.
     */
    private VCFHeaderLineCount getCountType(String key) {
        Object header = getHeader();
        if (header instanceof VCFHeader) {
            VCFInfoHeaderLine line = ((VCFHeader) header).getInfoHeaderLine(key);
            if (line != null) {
                return line.getCountType();
            }
        }
        return VCFHeaderLineCount.UNBOUNDED;
    }

    /**
     * Combine the elements of a value from {@code first} on.  Null if none of them is a number.  Elements that
     * are not numbers, such as the missing marker, are skipped rather than making the whole value unusable.
     */
    static Double aggregate(String[] parts, int first, Aggregation aggregation) {

        Double result = null;

        for (int i = first; i < parts.length; i++) {
            double d;
            try {
                d = Double.parseDouble(parts[i].trim());
            } catch (NumberFormatException e) {
                continue;
            }
            if (result == null) {
                result = d;
            } else if (aggregation == Aggregation.SUM) {
                result += d;
            } else {
                result = Math.max(result, d);
            }
        }
        return result;
    }

    /**
     * Normalize an attribute value for use as a color table key.  Multi-valued attributes are returned by htsjdk
     * as a list ("[a, b]"), which is reduced to "a,b" here so the key matches what the user sees.  Returns null
     * for a value that means "missing", which is drawn in the no-value color.
     */
    private static String normalizeAttributeValue(String value) {
        if (value == null) {
            return null;
        }
        String v = value.trim();
        if (v.startsWith("[") && v.endsWith("]")) {
            v = v.substring(1, v.length() - 1);
        }
        v = v.replaceAll("\\s*,\\s*", ",").trim();
        // "." is the VCF missing value marker.  htsjdk passes it through for a string attribute, and treated as
        // a value it would take a color of its own and appear in the legend as if it meant something.
        return v.isEmpty() || ".".equals(v) || "null".equals(v) ? null : v;
    }

    /**
     * Return the INFO header lines that make sense to color by, sorted by ID.  Numeric attributes are included:
     * selecting one asks for a color scale rather than coloring each distinct value (see
     * VariantTrackMenuHelper.defineScaleIfNeeded), so they must be offered.  Also included is any attribute a
     * color scheme covers, whatever its declared type.
     */
    public List<VCFInfoHeaderLine> getColorableInfoFields() {
        Object header = getHeader();
        if (!(header instanceof VCFHeader)) {
            return Collections.emptyList();
        }
        Set<String> schemeKeys = VariantColorSchemes.getKeys();
        return ((VCFHeader) header).getInfoHeaderLines().stream()
                .filter(line -> COLORABLE_INFO_TYPES.contains(line.getType())
                        || schemeKeys.contains(line.getID().toUpperCase()))
                .sorted(Comparator.comparing(VCFInfoHeaderLine::getID, String.CASE_INSENSITIVE_ORDER))
                .collect(Collectors.toList());
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
        // The variant band is not a "row", its height is independent of the display mode
        return VARIANT_BAND_HEIGHT;
    }

    public enum ColorMode {
        GENOTYPE, METHYLATION_RATE, ALLELE_FREQUENCY, NONE, ALLELE_FRACTION, ATTRIBUTE
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
            CircularViewUtilities.addVariants(svFeatures, getName(), CIRC_VIEW_DEFAULT_COLOR);
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

    /**
     * SQUISHED and EXPANDED set the genotype band (row) height.  The legacy COLLAPSED mode is translated to
     * "hide genotypes".
     *
     * @param mode
     */
    @Override
    public void setDisplayMode(DisplayMode mode) {
        if (mode == DisplayMode.COLLAPSED) {
            this.showGenotypes = false;
            mode = DisplayMode.EXPANDED;
        }
        super.setDisplayMode(mode);
    }

    /**
     * Legacy sessions store an explicit "squishedHeight".  Honor it if the track is in squished mode.
     */
    private void applyLegacySquishedHeight(int squishedHeight) {
        if (squishedHeight > 0 && getDisplayMode() == DisplayMode.SQUISHED) {
            setRowHeight(squishedHeight);
        }
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("showGenotypes")) {
            this.showGenotypes = Boolean.parseBoolean(element.getAttribute("showGenotypes"));
        }
        if (element.hasAttribute("squishedHeight")) {
            applyLegacySquishedHeight(Integer.parseInt(element.getAttribute("squishedHeight")));
        }
        if (element.hasAttribute("genotypeColorMode")) {
            this.genotypeColorMode = ColorMode.valueOf(element.getAttribute("genotypeColorMode"));
        } else if (element.hasAttribute("coloring")) {
            // backward compatibility
            this.genotypeColorMode = ColorMode.valueOf(element.getAttribute("coloring"));
        }

        if (element.hasAttribute("colorByAttribute")) {
            this.colorByAttribute = element.getAttribute("colorByAttribute");
            if (element.hasAttribute("attributeColorTable")) {
                getAttributeColorTable(colorByAttribute).restoreMapFromString(element.getAttribute("attributeColorTable"));
            }
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
        if (genotypeColorMode != ColorMode.GENOTYPE) {
            json.put("genotypeColorMode", genotypeColorMode.toString());
        }
        if (siteColorMode != null) {
            json.put("siteColorMode", siteColorMode.toString());
        }
        if (colorByAttribute != null) {
            json.put("colorByAttribute", colorByAttribute);

            // Colors the user chose.  Shaped like the igv.js "colorTable" track property.
            Map<String, Color> overrides = getAttributeColorOverrides(colorByAttribute);
            if (!overrides.isEmpty()) {
                org.json.JSONObject colorTable = new org.json.JSONObject();
                for (Map.Entry<String, Color> entry : overrides.entrySet()) {
                    colorTable.put(entry.getKey(), ColorUtilities.colorToString(entry.getValue()));
                }
                json.put("colorTable", colorTable);
            }

            // Colors IGV assigned from the palette.  Persisted so they are reproducible, they are otherwise
            // assigned in the order values are encountered.
            PaletteColorTable paletteColors = attributeColorTables.get(colorByAttribute);
            if (paletteColors != null && !paletteColors.getColorMap().isEmpty()) {
                json.put("attributeColorTable", paletteColors.getMapAsString());
            }
        }
    }


    @Override
    public void unmarshalJSON(org.json.JSONObject json) {

        super.unmarshalJSON(json);

        if (json.has("showGenotypes")) {
            this.showGenotypes = json.getBoolean("showGenotypes");
        }

        if (json.has("squishedHeight")) {
            applyLegacySquishedHeight(json.getInt("squishedHeight"));
        }

        if (json.has("genotypeColorMode")) {
            this.genotypeColorMode = ColorMode.valueOf(json.getString("genotypeColorMode"));
        }

        // Restore the attribute and its color table before the color mode, setColorByAttribute would override it
        if (json.has("colorByAttribute")) {
            this.colorByAttribute = json.getString("colorByAttribute");
            if (json.has("attributeColorTable")) {
                getAttributeColorTable(colorByAttribute).restoreMapFromString(json.getString("attributeColorTable"));
            }
            if (json.has("colorTable")) {
                org.json.JSONObject colorTable = json.getJSONObject("colorTable");
                for (String value : colorTable.keySet()) {
                    setAttributeColorOverride(colorByAttribute, value,
                            ColorUtilities.stringToColor(colorTable.getString(value)));
                }
            }
        }

        if (json.has("siteColorMode")) {
            this.siteColorMode = ColorMode.valueOf(json.getString("siteColorMode"));
        }
    }

}
