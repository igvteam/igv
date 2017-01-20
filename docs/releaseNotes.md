## **IGV 2.3,** released April 2013

**New Features and Improvements TESTING**

*   **Motif finder.** This new feature allows you to search for a particular nucleotide sequence in the reference genome, using either regular expression syntax or IUPAC ambiguity codes. Search results are displayed in two new tracks; one displays matches on the positive strand, and the other displays matches on the negative strand. Read more about it [here](http://www.broadinstitute.org/software/igv/node/251).

*   **Copy details of alignment coverage track to the clipboard. **If you hover the mouse cursor over a sequence alignment coverage track, the pop-up window shows details about the reads at that locus, including the number of reads, and the distribution of called nucleotides. You can now copy this information to your computer's clipboard by right-clicking and selecting _Copy Details to Clipboard _from the menu.

*   **Improved memory performance for BAM files.** In this release, sequence alignments that are no longer in view are purged from memory more frequently than before. "Out of memory" messages should happen less often.

*   **Gitools integration.** Gitools is a framework, developed at the Biomedical Research Park in Barcelona (PRBB), for analysis and visualization of genomic data. Data and results are represented as browsable heatmaps. Data can be exported from IGV in gitools format, or loaded directly into a running gitools session. Available under the _Tools_ menu. See http://www.gitools.org for details. 

##### Bug Fixes

*   Improved parsing of GFF3 files, so alternate coding sequences are displayed properly.
*   Fixed loading of BAM files over FTP.
*   Fixed bug where "Load Genome" dialog was sometimes blank.
*   Set autoscale to false when scale manually set, don't disable "Set Data Range" when autoscale set.
*   Trim any leading bases from VCF files which match the reference, so variants are displayed where the variation actually starts.

# **IGV 2.3.2 ** April 22, 2013

<div>

##### Bug Fixes

*   Cannot load large tdf or tabix files over http from some servers.

# **IGV 2.3.3 ** May 1, 2013

<div>

##### Bug Fixes

*   Minimum junction setting is ignored in splice junction track.
*   "Goto mate region" jumps to whole chromosome view when mate is on another chromosome.

# **IGV 2.3.4 ** May 13, 2013

##### Bug Fixes

*   File handes for "tribble" indexed files are not being closed.
*   Dragging in split screen / gene list view is broken.
*   Coverage track height cannot be changed.

# **IGV 2.3.6 ** May  2013

##### Bug Fixes

*   Null pointer exception creating snapshots from batch script.

</div>

</div>

# **IGV 2.3.7 ** May 30, 2013

<div>

*   The region-of-interest (ROI) batch command now takes an optional "name" argument.

##### Bug Fixes

*   When viewing VCF files with large numbers of samples view is not refreshed after scrolling vertically.
*   IGV fails on startup with Java 8

# **IGV 2.3.8 ** May 31, 2013

<div>

##### Bug Fixes

*   Fails to use proxy settings when performing byte range requests.
*   Searching for full chromosome causes no view change if already zoomed in

# **IGV 2.3.10  July 7**, 2013

<div>

##### **New features and Improvements**

*   Display message when BEDTools plugin encounters an error
*   Show genomic location in tooltip of alignment tracks
*   Sashimi plot for RNA-Seq alignments
*   Option to override user preferences via a batch comman>
*   Option to direct igvtools output to stdout for count and sort commands, using the special string 'stdout' as the output file.

##### **Bug Fixes**

*   igvtools outputs incorrect wig files (multiple track lines)
*   SVG files are not readable by some programs (e.g. gimp) on Windows platform due to incorrect encoding.

</div>

</div>

</div>

# **IGV 2.3.12  July 15**, 2013

<div>

##### **Bug Fixes**

*   Some bigwig files fail to load with a "null pointer exception"

# **IGV 2.3.13 August 2**, 2013

*   Pickup default  proxy settings from environment variables.

<div>**Bug Fixes**

*   UTRs not recognized in GFF files

# **IGV 2.3.14 August 6**, 2013

*   Implement autoscale for splice junction tracks

<div>**Bug Fixes**

*   Loading bedgraph files fail with "array index out of bounds" error

# **IGV 2.3.15 August 19**, 2013

<div>**Bug Fixes**

*   Show read sequence of unmapped mates 
*   Exons not associated with transcripts in GTF files
*   Splice junction tooltip incomplete in collapsed mode

</div>

</div>

</div>

</div>

# **IGV 2.3.16 August 20**, 2013

*   Add "export consensus sequence" option to alignment track menu

<div>

# **IGV 2.3.17 August 26**, 2013

<div>

*   Disable native Mac dialog when IGV is launched from Java Web Start 

# **IGV 2.3.18 September 5**, 2013

<div>

**Bug Fixes**

*   Autoscale option is reset after navigating to new locus
*   Reference sequence track is not shown from session files if there is no assoicated annotation track

# **IGV 2.3.19 September 27**, 2013

<div>

**Bug Fixes**

*   Popup messages appear when running batch scripts from the command line
*   Starting from the command line with a session file fails if no genome was previously selected
*   IGV hangs when viewing soft-clipped alignments at coordinates < 0  (5' end of contig)
*   Not all exons are displayed from Wormbase GFF files

# **IGV 2.3.23 October 24**, 2013

**New Features and Improvements**

*   Add preference to control display of default sample/track attributes
*   Performance improvements for loading BAM files over the web
*   Development version only: Remove tracks via delete or backspace key
*   Allow multiple motifs to be entered at once in MotifFinder
*   Allow grouping of VCF tracks

**Bug Fixes**

*   Starting igv with the -g option causes a null pointer exception when genome is in current directory
*   Cannot read BAM files when server doesn't return Content-Length
*   Cannot read eQTL files with missing q-value

# **IGV 2.3.26 January 17**, 2014

</div>

</div>

</div>

</div>

**New Features and Improvements**

*   Add option to [overlay data tracks](https://www.broadinstitute.org/software/igv/PopupMenus#DataTrack)
*   A new option for loading files from the ENCODE project has been added to the file menu for genome assemblies hg19, b37, and mm9.
*   Provide 1-click download of genome and sequence files. 
*   Add additional options for loading files by URL, can specify the index file and TDF file location with "index" and "coverage" query parameters respectively.
*   Allow EPS snapshot output, through use of EPSGraphics library. See [https://www.broadinstitute. <wbr> org/software/igv/third_party_ <wbr> tools#epsgraphics](https://www.broadinstitute.org/software/igv/third_party_tools#epsgraphics) for instructions
*   Add "remove" batch command, syntax "remove <trackName>"
*   Add 'ctrl-s' hotkey for sorting alignments by most recent sort
*   Read pairs are now downsampled together (both kept or both removed)

<wbr> <wbr>

**Bug Fixes**

*   Fallback from HEAD to GET for HTTP requests if HEAD fails. Useful for loading from cloud providers such as AWS
*   When long error messages are displayed, they are scrollable

<div>

# **IGV 2.3.46 March 16**, 2015

**New Features and Improvements**

*   **BLAT. **You can now do a BLAT search from a _feature_, _alignment_, or _region of interest, _of up to 8 kb in length as detailed [here](http://www.broadinstitute.org/software/igv/BLAT).
    *   **Features:** <span>Right-click on the feature in the track to see the BLAT command in the pop-up menu. The BLAT input sequence is the section of the reference genome defined by the feature start and end bounds.</span>
    *   **Alignments:** Right-click on the aligned read to see the BLAT command in the pop-up menu. The BLAT input sequence is the read sequence. It is _not _the sequence of the reference genome where the read was aligned.
    *   **Regions of Interest: **After [creating a region of interest](http://www.broadinstitute.org/software/igv/regionsofinterest), click on the region's red bar to see the BLAT command in the popup menu. The BLAT input sequence is the sequence of the reference genome defined by the region bounds.<span>   </span>

<span>​</span>By default, the BLAT server hosted at the UCSC Genome Browser is used, but this can be changed in the _Advanced_ preferences (_View > Preferences > Advanced_). Most UCSC derived genomes are supported, including human and mouse genomes.

*   <span>​</span>**Combine numeric tracks.** From the _Tools _menu, you can now dynamically create new tracks by combining two existing numeric tracks. Operators include add, subtract, multiply, and divide. Note that this does not overlay the two tracks, but rather creates a new data track that combines the data values from both tracks. For example, if you multiply two tracks, and one track has a data value of 10 at a particular locus, and the other track has a value of 2 at that locus, then the new track will have a value of 20. 

*   **Support for Supplementary Alignments**
    *   <span>Added a new alignment preference to filter supplementary alignments. See _Filter and shading options_ in </span>_View > Preferences > Alignments_<span>.</span>
    *   <span>Added _Group by supplementary flag_ to the alignment popup track menu.</span>
    *   <span>Enhanced the pop-up information panel, also referred to as a tooltip, for supplementary alignments.</span>

</div>

*   <span>**Misc**</span>
    *   <span>Session file paths are now relative where possible (i.e. data file path is relative to session path).</span>
    *   <span>Added a new batch command to go to whole genome view: “goto all”.  </span>
    *   <span>You can now optionally disable quality weighting of the allele fraction when displaying single nucleotide mismatches in the coverage track for sequencing data. See "Filter and shading options" in </span>_View > Preferences > Alignments_<span>.</span>
    *   <span>Nucleotide colors for the sequence track are now settable by adding or editing properties COLOR.A, COLOR.C, COLOR.T, COLOR.G, and COLOR.N in the IGV _pref.properties_ file (see the _igv _directory in your home directory). Similarly for a BAM track with properties SAM.COLOR.A, SAM.COLOR.C, etc., as outlined [here](http://www.broadinstitute.org/software/igv/prefs.properties).</span>
    *   <span>The exon numbering scheme has been changed to begin at first _coding_ exon, rather than first exon.</span>
    *   <span>Increased the keyboard left/right scroll speed to 50 pixels per key press.</span>
    *   <span>Improved memory management for smaller BAM file footprint.</span>

# **IGV 2.3.47 March 26**, 2015

**New Features and Improvements**

*   Add option to filter junctions on Sashimi plot by strand

**Bug Fixes**

*   Fix SSL certification errors when running under Java 1.8

# **IGV 2.3.49 April 10**, 2015

**Bug Fixes**

*   "Show centerline" ignored if alignment tracks are zoomed out

# **IGV 2.3.50 April 15**, 2015

**New Features and Improvements**

*   Option to add gene lists in "bed" format
*   Add  "count as pairs" and "extension factor" options to igvtools input form
*   Add option to scale fonts for very high resolution displays.  If the screen resolution is > 1.5 times the design resolution of 96 dpi all fonts are scaled proportionally.  Option is set on the "Genera" tab of the preferences window, default is off (don't scale).

**Bug Fixes**

*   Null pointer exception preventing startup when scaling fonts on some machines with very-high resolution screens.

# **IGV 2.3.51 April 20, 2015**

**Bug Fixes**

*   Null pointer exception preventing startup when scaling fonts on some machines with very-high resolution screens.

# **IGV 2.3.52 April 20**, 2015

**Bug Fixes**

*   BAM files fail to load and hang IGV if server returns a "302" status code for file not found.

# **IGV 2.3.53 May 26, 2015**

**Bug Fixes**

*   Bug -- #gffTags ignored if followed by track line
*   Concatenate bases of leading soft-clip operators

**New Features and Improvements**

*   Support for UCSC snp files (extension ".snp").
*   Check of active port listener before attempting OAuth.
*   Determine content length of remote file from partial content response header.
*   Print version update notice to log instead of popup dialog to prevent batch runs from failing.

# **IGV 2.3.55 May 28, 2015**

**Bug Fixes**

*   Bug -- ".bam.bai" index file not found in some situations

# **IGV 2.3.56 June 19 , 2015**

**New Features and Improvements**

*   Support for proxy "white list"
*   Improve treatment of cigar operators with padding
*   Support ".bdg" as a bedgraph file extension
*   Ignore "YC" tag in bam records if not convertible to a color

# **IGV 2.3.57 June 19 , 2015**

**Bug Fixes**

*   Check for "null" gap types when computing splice junctions.

# **IGV 2.3.58 August 4, 2015**

**New Features and Improvements**

*   Support for Encode gapped peak files (extension .gappedPeak).
*   Update htsjdk to version 1.138

# **IGV 2.3.59 August 7 , 2015**

**Bug Fixes**

*   Splice gaps in alignments sometimes drawn as deletion when show soft clips is on

**New Features and Improvements**

*   User preference for customizing VCF colors
*   User preference for turning on "show all bases"

# **IGV 2.3.60 September 21 , 2015**

**Bug Fixes**

*   Center position off-by-one when using "goto"

**New Features and Improvements**

*   New option to render bed files as arcs connecting the start and end position. To use specify "graphType=arc" on the track line. An initial height of at least 250 pixels is recommended (track line option height=250)
*   New bisulfite option to show all "C"s irrespective of context (context = None)
*   Session files now uses absolute resource paths by default. There is a new option to force use of relative paths in the "General" category of user preferences.

# **IGV 2.3.61 October 20 , 2015**

**Bug Fixes**

*   VCF and bed files not loadable from Google cloud storage (gs:// urls)

# **IGV 2.3.67 November 27 , 2015**

**Bug Fixes**

*   Recognize GFF version directives with major/minor version numbers

**New Features and Improvements**

*   Add support for exporting features from the Splice Junction track

# **IGV 2.3.68 January 13, 2016**

**New Features and Improvements**

*   VCF colors can now be customized. See the "Variants" tab of the Preferences dialog
*   New options for loading BAM files - alignments can now be hidden, showing coverage and/or splice junctions only. This option leads to significant memory savings when viewing individual alignments is not neccessary. See the new track display options on the "Alignments" tab of the Preferences dialog. These options control initial display upon loading a BAM file. Alignments can still be loaded on demand for specific loci by right-clicking on the coverage or splice-junction track and selecting _Show Alignment Track_. NOTE: When the Alignment track is off consider selecting "Display all tracks in a single panel" from the "General" tab of the Preferences dialog.
*   Add new menu item to autoscale multiple tracks as a group. Item is enabled if multiple numeric tracks are selected.

<a name="2.3.69"></a>

# **IGV 2.3.69 March 8, 2016**

**New Features and Improvements**

*   Support indexing of .mut and .maf mutation files with igvtools.
*   Add "Export Track Names" to right-click menu. Selected track names and associated meta-data (sample information) are exported to a file.
*   Support for loading Bionano .smap files.
*   **Note: Important change to VCF variant coloring.** VCF variants can now be colored by either allele fraction (AF) or allele fraction (computed from AC and AN fields). This option is specified on the Variant tab of the user preferences window. **_Calculations of Allele Fraction from genotype records is no longer performed in IGV_**, either AF or AC & AN must be specified in the INFO column to color variants.

**Bug Fixes**

*   Suppress confirm dialogs when running a batch script.
*   Fix problems with VCF allele counts.

<a name="2.3.71"></a>

# **IGV 2.3.71 March 21, 2016**

**Bug Fixes**

*   Cannot load a bam file with "refgene" in the name.
*   .mut mutation files fail to load with number format error.

<a name="2.3.72"></a>

# **IGV 2.3.72 April 1, 2016**

**Bug Fixes**

*   Alignments starting with deletions not rendered correctly.

<a name="2.3.73"></a>

# **IGV 2.3.73 May 25, 2016**

**New Features and Improvements**

*   igvtools now support sorting and indexing BAM files.
*   New MAF -> SAM converter available in igvtools.
*   Support group autoscale across panels.
*   Support group autoscale of alignment coverage tracks.
*   New "blat" item in Tools menu for blatting arbitrary sequence.

**Bug Fixes**

*   setMaxHeight command in batch files ignored on igv startup.

<a name="2.3.74"></a>

# **IGV 2.3.74 May 25, 2016**

**Bug Fixes**

*   Deletions starting out-of-view not rendered.

<a name="2.3.75"></a>

# **IGV 2.3.75 May 31, 2016**

**Bug Fixes**

*   Autoscale not working for coverage tracks.

<a name="2.3.76"></a>

# **IGV 2.3.76 June 9, 2016**

**New Features and Improvements**

*   Attributes from all constituitve parts of GFF features are now displayed in popup text.

**Bug Fixes**

*   Small chromosomes skipped when jumping through features (ctrl-f and ctrl-b).

<a name="2.3.78"></a>

# **IGV 2.3.78 June 28, 2016**

**New Features and Improvements**

*   Check sequence lengths in bam files for maximum supported length of 2^29-1.
*   Improved feature export for gff tracks.
*   Remove restriction on number of stacked features at a single locus.

<a name="2.3.79"></a>

# **IGV 2.3.79 July 9, 2016**

**New Features and Improvements**

*   Improved tool-tip text for alignment tracks.

**Bug Fixes**

*   OAuth credentials not passed for GA4GH readgroupsets.
*   Cannot specify index URL for google cloud storage BAM files.

<a name="2.3.80"></a>

# **IGV 2.3.80 August 4, 2016**

**New Features and Improvements**

*   Support Google gs:// urls from batch scripts.
*   Enable loading ga4gh read group sets from batch scripts.
*   Support for setting igv directory from command line switch (--igvDirectory).
*   Display error dialog if fasta index file cannot be created.
*   Allow setting oAuth access token via a port command.

**Bug Fixes**

*   Copy number segmenets not always drawn if start is close to left edge of screen.

<a name="2.3.81"></a>

# **IGV 2.3.81 August 29, 2016**

**New Features and Improvements**

*   New "color-by library" option for alignment tracks
*   Session loading by drag&drop

**Bug Fixes**

*   Relative genome paths in session files don't work.

<a name="2.3.82"></a>

# **IGV 2.3.82 September 27, 2016**

**Bug Fixes**

*   Sashimi plots do not autoscale.

<a name="2.3.83"></a>

# **IGV 2.3.83 September 27, 2016**

**Bug Fixes**

*   Create new panel for bam files loaded from batch script (git issue #322)
*   Fix git issue #317\. Server not loading data from tracks where the index file has a distinct URL

<a name="2.3.84"></a>

# **IGV 2.3.84 October 27, 2016**

**Bug Fixes**

*   Fix: GenomeSpace paths with spaces will not load.

<a name="2.3.85"></a>

# **IGV 2.3.85 October 30, 2016**

**Bug Fixes**

*   Line plots artificats in whole genome view from TDF files

**New Features and Improvements**

*   Support "points" plot type (with window function == none) for whole genome view from TDF files

<a name="2.3.86"></a>

# **IGV 2.3.86 November 1, 2016**

**Bug Fixes**

*   Incorrect (black & white) colors for "seg" copy number files.
