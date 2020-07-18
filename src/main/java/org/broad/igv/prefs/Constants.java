package org.broad.igv.prefs;

import java.util.Arrays;

/**
 * Created by jrobinso on 1/20/17.
 */
final public class Constants {

    private Constants() {}  // Prevent instantiation

    // Preference sets
    public static final String NULL_CATEGORY = "NULL";
    public static final String THIRD_GEN = "THIRD_GEN";
    public static final String RNA = "RNA";

    // UI behavior
    public static final String TOOLTIP_INITIAL_DELAY = "TOOLTIP.INITIAL_DELAY";
    public static final String TOOLTIP_RESHOW_DELAY = "TOOLTIP.RESHOW_DELAY";
    public static final String TOOLTIP_DISMISS_DELAY = "TOOLTIP.DISMISS_DELAY";
    public static final String DEFAULT_FONT_SIZE = "DEFAULT_FONT_SIZE";
    public static final String DEFAULT_FONT_FAMILY = "DEFAULT_FONT_FAMILY";
    public static final String DEFAULT_FONT_ATTRIBUTE = "DEFAULT_FONT_ATTRIBUTE";
    public static final String ENABLE_ANTIALISING = "ENABLE_ANTIALIASING";
    public static final String SCALE_FONTS = "SCALE_FONTS";
    public static final String NAME_PANEL_WIDTH = "NAME_PANEL_WIDTH";
    public static final String BACKGROUND_COLOR = "BACKGROUND_COLOR";
    public static final String SHOW_ATTRIBUTE_VIEWS_KEY = "IGV.track.show.attribute.views";
    public static final String SHOW_SINGLE_TRACK_PANE_KEY = "IGV.single.track.pane";
    public static final String DETAILS_BEHAVIOR_KEY = "DETAILS_BEHAVIOR";
    public static final String SKIP_VERSION = "SKIP_VERSION";
    public static final String SHOW_DEFAULT_TRACK_ATTRIBUTES = "SHOW_DEFAULT_TRACK_ATTRIBUTES";
    public static final String SHOW_GENOME_SERVER_WARNING = "SHOW_GENOME_SERVER_WARNING";
    public static final String CN_FREQ_AMP_THRESHOLD = "CN_FREQ.AMP_THRESHOLD";
    public static final String CN_FREQ_DEL_THRESHOLD = "CN_FREQ.DEL_THRESHOLD";

    //
    public static final String RECENT_SESSIONS = "IGV.Session.recent.sessions";
    public static final String LAST_EXPORTED_REGION_DIRECTORY = "LAST_EXPORTED_REGION_DIRECTORY";
    public static final String LAST_TRACK_DIRECTORY = "LAST_TRACK_DIRECTORY";
    public static final String LAST_SNAPSHOT_DIRECTORY = "LAST_SNAPSHOT_DIRECTORY";
    public static final String LAST_GENOME_IMPORT_DIRECTORY = "LAST_GENOME_IMPORT_DIRECTORY";
    public static final String DEFINE_GENOME_INPUT_DIRECTORY_KEY = "DEFINE_GENOME_INPUT_DIRECTORY_KEY";
    public static final String DEFAULT_GENOME = "DEFAULT_GENOME_KEY";
    public static final String AUTO_UPDATE_GENOMES = "AUTO_UPDATE_GENOMES";
    public static final String FRAME_BOUNDS_KEY = "IGV.Bounds";


    public static final String GENOMES_SERVER_URL = "IGV.genome.sequence.dir";
    public static final String BLAT_URL = "BLAT_URL";
    public static final String EXTVIEW_URL = "EXTVIEW_URL";
    public static final String DATA_SERVER_URL_KEY = "MASTER_RESOURCE_FILE_KEY";

    public static final String CRAM_CACHE_SEQUENCES = "CRAM.CACHE_SEQUENCES";
    public static final String CRAM_CACHE_DIRECTORY = "CRAM.CACHE_DIRECTORY";
    public static final String CRAM_CACHE_SIZE = "CRAM.CACHE_SIZE";

    // Search ("go to") options
    public static final String SEARCH_ZOOM = "SEARCH_ZOOM";
    public static final String FLANKING_REGION = "FLANKING_REGION";

    // Generic track options
    public static final String BYPASS_FILE_AUTO_DISCOVERY = "BYPASS_FILE_AUTO_DISCOVERY";
    public static final String TRACK_ATTRIBUTE_NAME_KEY = "TRACK_ATTRIBUTE_NAME_KEY";
    public static final String INITIAL_TRACK_HEIGHT = "15";
    public static final String COLOR_SCALE_KEY = "COLOR_SCALE_";
    public static final String TRACK_HEIGHT_KEY = "IGV.track.height";
    public static final String CHART_TRACK_HEIGHT_KEY = "IGV.chart.track.height";
    public static final String CHART_SHOW_ALL_HEATMAP = "CHART.SHOW_ALL_HEATMAP";
    public static final String SHOW_REGION_BARS = "SHOW_REGION_BARS";
    public static final String DEFAULT_VISIBILITY_WINDOW = "DEFAULT_VISIBILITY_WINDOW";
    public static final String EXPAND_FEAUTRE_TRACKS = "EXPAND_FEATURE_TRACKS";
    public static final String IGV_PLUGIN_LIST_KEY = "IGV_PLUGIN_LIST";
    public static final String SASHIMI_SHOW_COVERAGE = "SASHIMI.SHOW_COVERAGE";
    public static final String GENE_LIST_BED_FORMAT = "GENE_LIST_BED_FORMAT";
    public static final String SESSION_RELATIVE_PATH = "SESSION.RELATIVE_PATH";
    public static final String SHOW_SIZE_WARNING = "SHOW_SIZE_WARNING";

    // Interaction track
    public static final String ARC_TYPE = "ARC_TYPE";
    public static final String ARC_DIRECTION = "ARC_DIRECTION";
    public static final String ARC_BLOCKS = "ARC_BLOCKS";

    // Chart (bar, heatmap, plots) options
    public static final String NORMALIZE_COVERAGE = "NORMALIZE_COVERAGE";
    public static final String CHART_DRAW_TOP_BORDER = "CHART.DRAW_TOP_BORDER";
    public static final String CHART_DRAW_BOTTOM_BORDER = "CHART.DRAW_BOTTOM_BORDER";
    public static final String CHART_COLOR_BORDERS = "CHART.COLOR_BORDERS";
    public static final String CHART_DRAW_Y_AXIS = "CHART.DRAW_AXIS";
    public static final String CHART_DRAW_TRACK_NAME = "CHART.DRAW_TRACK_NAME";
    public static final String CHART_COLOR_TRACK_NAME = "CHART.COLOR_TRACK_NAME";
    public static final String CHART_AUTOSCALE = "CHART.AUTOSCALE";
    public static final String CHART_SHOW_DATA_RANGE = "CHART.SHOW_DATA_RANGE";

    // Alignment options
    public static final String SAM_ALLELE_THRESHOLD = "SAM.ALLELE_THRESHOLD";
    public static final String SAM_ALLELE_USE_QUALITY = "SAM.ALLELE_USE_QUALITY";
    public static final String SAM_QUALITY_THRESHOLD = "SAM.QUALITY_THRESHOLD";
    public static final String SAM_ALIGNMENT_SCORE_THRESHOLD = "SAM.ALIGNMENT_SCORE_THRESHOLD";
    public static final String SAM_COMPUTE_ISIZES = "SAM.COMPUTE_ISIZES";
    public static final String SAM_MAX_INSERT_SIZE_THRESHOLD = "SAM.INSERT_SIZE_THRESHOLD";
    public static final String SAM_MIN_INSERT_SIZE_THRESHOLD = "SAM.MIN_INSERT_SIZE_THRESHOLD";
    public static final String SAM_MAX_INSERT_SIZE_PERCENTILE = "SAM.ISIZE_MAX_PERCENTILE";
    public static final String SAM_MIN_INSERT_SIZE_PERCENTILE = "SAM.MIN_ISIZE_MIN_PERCENTILE";
    public static final String SAM_AUTO_SORT = "SAM.AUTOSORT";
    public static final String SAM_SHADE_CENTER = "SAM.SHADE_CENTER";
    public static final String SAM_SHOW_CENTER_LINE = "SAM.SHOW_CENTER_LINE";
    public static final String SAM_SHOW_REF_SEQ = "SAM.SHOW_REF_SEQ";
    public static final String SAM_SHOW_COV_TRACK = "SAM.SHOW_COV_TRACK";
    public static final String SAM_SHADE_BASES = "SAM.SHADE_BASE_QUALITY";
    public static final String SAM_BASE_QUALITY_MIN = "SAM.BASE_QUALITY_MIN";
    public static final String SAM_BASE_QUALITY_MAX = "SAM.BASE_QUALITY_MAX";
    public static final String SAM_FILTER_ALIGNMENTS = "SAM.FILTER_ALIGNMENTS";
    public static final String SAM_FILTER_SECONDARY_ALIGNMENTS = "SAM.FILTER_SECONDARY_ALIGNMENTS";
    public static final String SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS = "SAM.FILTER_SUPPLEMENTARY_ALIGNMENTS";
    public static final String SAM_FILTER_URL = "SAM.FILTER_URL";
    public static final String SAM_HIDDEN_TAGS = "SAM.HIDDEN_TAGS";
    public static final String SAM_MAX_VISIBLE_RANGE = "SAM.MAX_VISIBLE_RANGE";
    public static final String SAM_FILTER_DUPLICATES = "SAM.FILTER_DUPLICATES";
    public static final String SAM_SHOW_DUPLICATES = "SAM.SHOW_DUPLICATES";
    public static final String SAM_QUICK_CONSENSUS_MODE = "SAM.QUICK_CONSENSUS_MODE";
    public static final String SAM_SHOW_SOFT_CLIPPED = "SAM.SHOW_SOFT_CLIPPED";
    public static final String SAM_MAX_SOFT_CLIP = "SAM.MAX_SOFT_CLIP";
    public static final String SAM_FLAG_UNMAPPED_PAIR = "SAM.FLAG_UNMAPPED_PAIR";
    public static final String SAM_SAMPLING_COUNT = "SAM.MAX_LEVELS"; // Sampling count
    public static final String SAM_SAMPLING_WINDOW = "SAM.SAMPLING_WINDOW";
    public static final String SAM_DOWNSAMPLE_READS = "SAM.DOWNSAMPLE_READS";
    public static final String SAM_SORT_OPTION = "SAM.SORT_OPTION";
    public static final String SAM_GROUP_OPTION = "SAM.GROUP_OPTION";
    public static final String SAM_SHOW_ALL_BASES = "SAM.SHOW_ALL_BASES";
    public static final String SAM_SHOW_MISMATCHES = "SAM.SHOW_MISMATCHES";
    public static final String SAM_COLOR_BY = "SAM.COLOR_BY";
    public static final String SAM_COLOR_BY_TAG = "SAM.COLOR_BY_TAG";
    public static final String SAM_SORT_BY_TAG = "SAM.SORT_BY_TAG";
    public static final String SAM_GROUP_BY_TAG = "SAM.GROUP_BY_TAG";
    public static final String SAM_LINK_BY_TAGS = "SAM.LINK_BY_TAGS";
    public static final String SAM_GROUP_BY_POS = "SAM.GROUP_BY_POS";
    public static final String SAM_BISULFITE_CONTEXT = "SAM.BISULFITE_CONTEXT";
    public static final String SAM_FILTER_FAILED_READS = "SAM.FILTER_FAILED_READS";
    public static final String SAM_FLAG_ZERO_QUALITY = "SAM.FLAG_ZERO_QUALITY";
    public static final String SAM_SHOW_JUNCTION_TRACK = "SAM.SHOW_JUNCTION_TRACK";
    public static final String SAM_JUNCTION_MIN_FLANKING_WIDTH = "SAM.JUNCTION_MIN_FLANKING_WIDTH";
    public static final String SAM_JUNCTION_MIN_COVERAGE = "SAM.JUNCTION_MIN_COVERAGE";
    public static final String SAM_SHOW_INSERTION_MARKERS = "SAM.SHOW_INSERTION_MARKERS";
    public static final String SAM_SHOW_JUNCTION_FLANKINGREGIONS = "SAM.SHOW_JUNCTION_FLANKINGREGIONS";
    public static final String SAM_NOMESEQ_ENABLED = "SAM.NOMESEQ_ENABLED";
    public static final String SAM_COUNT_DELETED_BASES_COVERED = "SAM.COUNT_DELETED_BASES_COVERED";
    public static final String SAM_FLAG_LARGE_INDELS = "SAM.FLAG_LARGE_INDELS";
    public static final String SAM_LARGE_INDELS_THRESHOLD = "SAM.LARGE_INSERTIONS_THRESOLD";
    public static final String SAM_FLAG_CLIPPING = "SAM.FLAG_CLIPPING";
    public static final String SAM_CLIPPING_THRESHOLD = "SAM.CLIPPING_THRESHOLD";
    public static final String SAM_SHOW_GROUP_SEPARATOR = "SAM.SHOW_GROUP_SEPARATOR";
    public static final String SAM_REDUCED_MEMORY_MODE = "SAM.REDUCED_MEMORY_MODE";
    public static final String SAM_HIDE_SMALL_INDEL = "SAM.HIDE_SMALL_INDEL";
    public static final String SAM_SMALL_INDEL_BP_THRESHOLD = "SAM.SMALL_INDEL_BP_THRESHOLD";
    public static final String SAM_SHOW_ALIGNMENT_TRACK = "SAM.SHOW_ALIGNMENT_TRACK";
    public static final String SAM_COLOR_A = "SAM.COLOR.A";
    public static final String SAM_COLOR_C = "SAM.COLOR.C";
    public static final String SAM_COLOR_T = "SAM.COLOR.T";
    public static final String SAM_COLOR_G = "SAM.COLOR.G";
    public static final String SAM_COLOR_N = "SAM.COLOR.N";
    public static final String SAM_DISPLAY_MODE = "SAM.DISPLAY_MODE";
    public static final String KNOWN_SNPS = "KNOWN_SNPS_FILE";

    // Sequence track settings
    public static final String SHOW_SEQUENCE_TRANSLATION = "SHOW_SEQUENCE_TRANSLATION";
    public static final String MAX_SEQUENCE_RESOLUTION = "MAX_SEQUENCE_RESOLUTION";
    public static final String COLOR_A = "COLOR.A";
    public static final String COLOR_C = "COLOR.C";
    public static final String COLOR_T = "COLOR.T";
    public static final String COLOR_G = "COLOR.G";
    public static final String COLOR_N = "COLOR.N";

    // Variant (VCF) track settings
    public static final String VARIANT_COLOR_BY_ALLELE_FREQ = "VARIANT_COLOR_BY_ALLELE_FREQ";
    public static final String HOMREF_COLOR = "HOMREF.COLOR";
    public static final String HETVAR_COLOR = "HETVAR.COLOR";
    public static final String HOMVAR_COLOR = "HOMVAR.COLOR";
    public static final String NOCALL_COLOR = "NOCALL.COLOR";
    public static final String AF_REF_COLOR = "AF_REF.COLOR";
    public static final String AF_VAR_COLOR = "AF_VAR.COLOR";

    // Heatmap settings
    public static final String NO_DATA_COLOR = "NO_DATA.COLOR";
    public static final String NO_CALL_COLOR = "NO_CALL.COLOR";

    // "Mut" and "MAF" mutation track settings
    public static final String MUTATION_COLOR_TABLE = "MUTATION_COLOR_TABLE";
    public static final String MUTATION_INDEL_COLOR_KEY = "MUTATION_INDEL_COLOR_KEY";
    public static final String MUTATION_MISSENSE_COLOR_KEY = "MUTATION_MISSENSE_COLOR_KEY";
    public static final String MUTATION_NONSENSE_COLOR_KEY = "MUTATION_NONSENSE_COLOR_KEY";
    public static final String MUTATION_SPLICE_SITE_COLOR_KEY = "MUTATION_SPLICE_SITE_COLOR_KEY";
    public static final String MUTATION_SYNONYMOUS_COLOR_KEY = "MUTATION_SYNONYMOUS_COLOR_KEY";
    public static final String MUTATION_TARGETED_REGION_COLOR_KEY = "MUTATION_TARGETED_REGION_COLOR_KEY";
    public static final String MUTATION_UNKNOWN_COLOR_KEY = "MUTATION_UNKNOWN_COLOR_KEY";
    public static final String OVERLAY_MUTATION_TRACKS = "OVERLAY_TRACKS_KEY";
    public static final String SHOW_ORPHANED_MUTATIONS = "SHOW_ORPHANED_MUTATIONS";
    public static final String OVERLAY_ATTRIBUTE_KEY = "OVERLAY_ATTRIBUTE_KEY";
    public static final String OVERLAY_MUTATIONS_WHOLE_GENOME = "OVERLAY_MUTATIONS_WHOLE_GENOME";
    public static final String COLOR_MUTATIONS = "COVER_OVERLAY_KEY";

    // GWAS track options
    public static final String GWAS_TRACK_HEIGHT = "GWAS_TRACK_HEIGHT";
    public static final String GWAS_DESCRIPTION_CACHE_SIZE = "GWAS_DESCRIPTION_CACHE_SIZE";
    public static final String GWAS_MIN_POINT_SIZE = "GWAS_MIN_POINT_SIZE";
    public static final String GWAS_MAX_POINT_SIZE = "GWAS_MAX_POINT_SIZE";
    public static final String GWAS_USE_CHR_COLORS = "GWAS_USE_CHR_COLORS";
    public static final String GWAS_SINGLE_COLOR = "GWAS_SINGLE_COLOR";
    public static final String GWAS_ALTERNATING_COLORS = "GWAS_ALTERNATING_COLORS";
    public static final String GWAS_PRIMARY_COLOR = "GWAS_PRIMARY_COLOR";
    public static final String GWAS_SECONDARY_COLOR = "GWAS_SECONDARY_COLOR";
    public static final String GWAS_SHOW_AXIS = "GWAS_SHOW_AXIS";

    // Gene expression track options
    public static final String PROBE_MAPPING_KEY = "PROBE_MAPPING_KEY";
    public static final String PROBE_MAPPING_FILE = "PROBE_MAPPING_FILE";
    public static final String USE_PROBE_MAPPING_FILE = "USE_PROBE_MAPPING_FILE";

    // Genome space
    public static final String GENOME_SPACE_ENABLE = "GENOME_SPACE_ENABLE";
    public static final String GENOME_SPACE_DM_SERVER = "GENOME_SPACE_DM_SERVER";
    public static final String GENOME_SPACE_ATM_SERVER = "GENOME_SPACE_ATM_SERVER";
    public static final String GENOME_SPACE_IDENTITY_SERVER = "GENOME_SPACE_IDENTITY_SERVER";

    // Google
    public static final String GOOGLE_API_KEY = "GOOGLE_API_KEY";
    public static final String GOOGLE_PROJECT = "GOOGLE_PROJECT";
    public static final String ENABLE_GOOGLE_MENU = "ENABLE_GOOGLE_MENU";
    public static final String SAVE_GOOGLE_CREDENTIALS = "SAVE_GOOGLE_CREDENTIALS";

    // CBIO connections
    public static final String CBIO_MUTATION_THRESHOLD = "CBIO_MUTATION_THRESHOLD";
    public static final String CBIO_AMPLIFICATION_THRESHOLD = "CBIO_AMPLIFICATION_THRESHOLD";
    public static final String CBIO_DELETION_THRESHOLD = "CBIO_DELETION_THRESHOLD";
    public static final String CBIO_EXPRESSION_UP_THRESHOLD = "CBIO_EXPRESSION_UP_THRESHOLD";
    public static final String CBIO_EXPRESSION_DOWN_THRESHOLD = "CBIO_EXPRESSION_DOWN_THRESHOLD";

    // Proxy settings
    public static final String USE_PROXY = "PROXY.USE";
    public static final String PROXY_HOST = "PROXY.HOST";
    public static final String PROXY_PORT = "PROXY.PORT";
    public static final String PROXY_AUTHENTICATE = "PROXY.AUTHENTICATE";
    public static final String PROXY_USER = "PROXY.USERNAME";
    public static final String PROXY_PW = "PROXY.PW";
    public static final String PROXY_TYPE = "PROXY.TYPE";
    public static final String PROXY_WHITELIST = "PROXY.WHITELIST";

    // Port settings
    public static final String PORT_ENABLED = "PORT_ENABLED";
    public static final String PORT_NUMBER = "PORT_NUMBER";

    // Database support -- never deployed
    public static final String DB_ENABLED = "DB_ENABLED";
    public static final String DB_HOST = "DB_HOST";
    public static final String DB_NAME = "DB_NAME";
    public static final String DB_PORT = "DB_PORT";

    // OAuth provisioning
    public static final String PROVISIONING_URL = "PROVISIONING.URL";

    // Experimental
    public static final String SCORE_VARIANTS = "SCORE_VARIANTS";




    /**
     * List of keys that affect the alignments loaded.  This list is used to trigger a reload, if required.
     * Not all alignment preferences need trigger a reload, this is a subset.
     */
    static java.util.List<String> SAM_RELOAD_KEYS = Arrays.asList(
            SAM_QUALITY_THRESHOLD,
            SAM_ALIGNMENT_SCORE_THRESHOLD,
            SAM_FILTER_ALIGNMENTS,
            SAM_FILTER_URL,
            SAM_MAX_VISIBLE_RANGE,
            SAM_SHOW_DUPLICATES,
            SAM_SHOW_SOFT_CLIPPED,
            SAM_SAMPLING_COUNT,
            SAM_SAMPLING_WINDOW,
            SAM_FILTER_FAILED_READS,
            SAM_DOWNSAMPLE_READS,
            SAM_FILTER_SECONDARY_ALIGNMENTS,
            SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS,
            SAM_JUNCTION_MIN_FLANKING_WIDTH,
            SAM_JUNCTION_MIN_COVERAGE
    );

    /**
     * List of keys that do not affect the alignments loaded but do affect how those
     * alignments are drawn.  A refresh is softer than a reload.
     */
    static java.util.List<String> SAM_REFRESH_KEYS = Arrays.asList(
            SAM_QUICK_CONSENSUS_MODE,
            SAM_ALLELE_THRESHOLD,
            SAM_FLAG_LARGE_INDELS,
            SAM_LARGE_INDELS_THRESHOLD,
            SAM_SHOW_INSERTION_MARKERS
    );

}
