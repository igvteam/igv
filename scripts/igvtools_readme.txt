This package contains command line utilities for preprocessing, computing
feature count density (coverage),  sorting, and indexing data files.
See also http://www.broadinstitute.org/software/igv/igvtools_commandline.

***************************************************************************
Java 11 is required for this release.  See our website for more
information about support for Java 8.
***************************************************************************

---------------------------------------------------------------------------
Starting with shell scripts
---------------------------------------------------------------------------
The utilities are invoked from one of the following scripts:

   igvtools (command line version for Linux and  macOS 10.x)
   igvtools_gui (gui version for Linux and  macOS 10.x)
   igvtools_gui_hidpi (HiDPI gui version for Linux)
   igvtools_gui.command (alternative double-clickable gui version for macOS 10.x)

   igvtools.bat (command line version for Windows)
   igvtools_gui.bat (gui version for Windows)

The general form of the command-line version is:

   igvtools [command] [options][arguments]
or
   igvtools.bat [command] [options][arguments]

Recognized commands, options,arguments, and file types are described below.

---------------------------------------------------------------------------
Starting with java
---------------------------------------------------------------------------

Igvtools can also be started directly using java as shown below.  This option
allows more control over java parameters, such as the maximum memory to
allocate.  In the example below igvtools is started with 1500 MB of memory
allocated and launched in the location where you have unpacked IGVTools.

   java -Xmx1500m --module-path=lib @igv.args --module=org.igv/org.broad.igv.tools.IgvTools [command] [options][arguments]

To start with a gui the command is

   java -Xmx1500m --module-path=lib @igv.args --module=org.igv/org.broad.igv.tools.IgvTools gui

Note that the command line has become more complex with Java 11 compared to
Java 8.  We recommend the shell scripts above for most users.

---------------------------------------------------------------------------
Memory settings
---------------------------------------------------------------------------

The scripts above allocate a fixed amount of memory.  If this amount is not
available on your platform you will get an obscure error along the lines of
"Could not start the Virtual Machine".   If this happens you will need to
reduce the amount of memory requested.  The memory is set via a "-Xmx" parameter. 
For example -Xmx1500m  requests 1500 MB,  -Xmx1g requests 1 gigabyte.


---------------------------------------------------------------------------
Other Java settings
---------------------------------------------------------------------------

You can override this memory setting (and other Java-related defaults) by editing 
IGV's java_arguments file, found here (create it if it doesn't exist):
   $HOME/.igv/java_arguments           (Mac and Linux)
   %USERPROFILE%/.igv/java_arguments   (Windows)

Specifically set the value of the "-Xmx" parameter, which will be commented-out (with
a '#' character) by default when this file is created.  For example, to start IGV with 
8 GB of memory, uncomment this line and set the value to 

   -Xmx8g

This will override the default 4GB memory specification.

Other Java-related command-line options can also be set in this file, though changing anything
beyond the memory specification is for advanced users only and is not recommended.  See
   https://docs.oracle.com/en/java/javase/11/tools/java.html
for more information on the Java 11 command line, and
   https://docs.oracle.com/en/java/javase/11/tools/java.html#GUID-4856361B-8BFD-4964-AE84-121F5F6CF111
in particular for specifics of the "java_arguments" file format.

---------------------------------------------------------------------------
HiDPI settings
---------------------------------------------------------------------------

HiDPI is supported natively by Java on Mac and Windows.  Users on these
platforms can ignore this section. 

For Linux users, the igvtools_gui_hidpi script is set up for 2x scaling.  
To modify it to do 4x scaling, for example, change the value

   -Dsun.java2d.uiScale=2

to

   -Dsun.java2d.uiScale=4

Fractional values are *NOT* supported at this time.  Note that here again, you can add
this specification to the java_arguments file instead of editing the launcher scripts.
Doing so will allow the standard 'igvtools_gui' to work properly on HiDPI screens without
the need for the 'igvtools_gui_hidpi' script.

---------------------------------------------------------------------------
Genome
---------------------------------------------------------------------------

The genome argument in the tile and count command can be either an id, or
a full path to an IGV .genome file.  The id for IGV supplied genomes are
listed below.  Genome definitions corresponding to these files are in the
"genomes" subdirectory of the igvtools install.  The id is derived by removing
the .extension from the filename.

---------------------------------------------------------------------------
COMMANDS
---------------------------------------------------------------------------

The recognized commands are tile, count, sort, and index.  Note that these
utilities are for working with ascii file formats, including SAM, but
do not work with BAM files.  For manipulating BAM files use samtools (http://samtools.sourceforge.net/).

---------------------------------------------------------------------------
Command "tile"
---------------------------------------------------------------------------
Warning: This command is deprecated. Use "toTDF" instead.

---------------------------------------------------------------------------
Command "toTDF"
---------------------------------------------------------------------------

The "toTDF" command converts a sorted data input file to  a binary tiled
data (.tdf) file. Input file formats supported  are .wig, .cn, .igv,
and .gct, TCGA mage-tab files, and "list" files.

List files are text files containing a list of files in one of the supported formats,
one file per line. When using a list file the format of the contained files must be
specified explicitly with the "fileType" parameter.  List files must end with the
extension ".list".  File paths can be absolute or relative to the directory containing
the list file.

Usage:

  igvtools toTDF [options]  [inputFile] [outputFile] [genome]


Required arguments:

  inputFile    The input file (see supported formats above).

  outputFile   Binary output file.  Must end in ".tdf".

  genome       A genome id or filename. See details below. Default is hg18.

Options:

  -z, --maxZoom num       Specifies the maximum zoom level to precompute. The default
               value is 7 and is sufficient for most files. To reduce file
               size at the expense of IGV performance this value can be
               reduced.

  -f, --windowFunctions  list     A comma delimited list specifying window functions to use
               when reducing the data to precomputed tiles.   Allowed
               values are  min, max,  mean, median, p2, p10, p90, and p98.
               The "p" values represent percentile, so p2=2nd percentile,
               etc.

  -p, --probeFile file      Specifies a "bed" file to be used to map probe identifiers
               to locations.  This option is useful when preprocessing gct
               files.  The bed file should contain 4 columns:
                  chr start end name
               where name is the probe name in the gct file.

  --fileType   Explicitly specify the file type.  This is a required parameter  for TCGA mage-tab and ".list" files.
               Possible values are mage-tab, .wig, .cn, .igv, and .gct.   Only mage-tab files downloaded from the
               TCGA data center or related sights are supported at this time.


  Conversion of ".gct" and "mage-tab" files results in the creation of an ".igv" file, which is sorted by genome
  position using the "sort" command.  For this case the following optional parameters can be specified.

  -t, --tmpDir tmpdir  Specify a temporary working directory.  For large input files
               this directory will be used to store intermediate results of
               the sort. The default is the users temp directory.

  -m, --maxRecords number  The maximum number of records to keep in memory during the
               sort.  The default value is 500000.  Increase this number
               if you receive "too many open files" errors.   Decrease it
               if you experience "out of memory" errors.



Example:

      igvtools toTDF -z 5  copyNumberFile.cn copyNumberFile.tdf hg18


Notes:

Data file formats, with the exception of .gct files, must be sorted by
start position.  If necessary files can be sorted with the "sort" command
described below.  Attempting to preprocess an unsorted file will result
in an  error.

---------------------------------------------------------------------------
Command "count"
---------------------------------------------------------------------------

The "count" command computes average feature density over a specified
window size across the genome. Common usages include computing coverage
for alignment files and counting hits in Chip-seq experiments. Supported
file formats are .sam,  .bam,  .aligned,  .sorted.txt,  and .bed, and
.bam.list files.  The latter format is a plain text file containing a list
of alignment or bed files, one file per line.

Usage:

  igvtools count [options] [inputFile] [outputFile] [genome]

Required arguments:

  inputFile    The input file (see supported formats above).

  outputFile   Either a binary tdf file, a text wig file, or both.  The output file type is determined
               by file extension, for example "output.tdf".  To output both formats supply two file names
               separated by a commas,  for example  "outputBinary.tdf,outputText.wig". To display feature
               intensity in IGV, the density must be computed with this command, and the resulting file
               must be named <feature track filename>.tdf.
               The special string "stdout" can be used in either position, in which case the output will
               be written to the standard output stream in wig format.

  genome       A genome id or filename. See details below. Default is hg18.

Options:

  -z, --maxZoom num       Specifies the maximum zoom level to precompute.

  -w, --windowSize num       The window size over which coverage is averaged. Defaults
               to 25 bp.

  -e, --extFactor num   The read or feature is extended by the specified distance
               in bp prior to counting. This option is useful for chip-seq
               and rna-seq applications. The value is generally set to the
               average fragment length of the library minus the average read length.


   --preExtFactor  num   The read is extended upstream from the 5' end by the specified distance.

   --postExtFactor num   Effectively overrides the read length, defines the downstream extent
               from the 5' end.  Intended for use with preExtFactor.


  -f, --windowFunctions  list     A comma delimited list specifying window functions to use
               when reducing the data to precomputed tiles.   Possible
               values are  min, max,  mean, median, p2, p10, p90, and p98.
               The "p" values represent percentile, so p2=2nd percentile,
               etc.

  --strands [arg] By default, counting is combined among both strands.
                This setting outputs the count for each strand separately.
                Legal argument values are 'read' or 'first'.
                'read' Separates count by 'read' strand, 'first' uses the first in pair strand.
                Results are saved in a separate column for .wig output, and a separate track
                for TDF output.

  --bases       Count the occurrence of each base (A,G,C,T,N). Takes no arguments.
                Results are saved in a separate column for .wig output, and a separate track for TDF output.
  
  --query [querystring] Only count a specific region. Query string has syntax <chr>:<start>-<end>. e.g. chr1:100-1000.
                        Input file must be indexed.
  
  --minMapQuality [mqual]   Set the minimum mapping quality of reads to include. Default is 0.

  --includeDuplicates    Include duplicate alignments in count. Default false.  If this flag is included, duplicates
                          are counted. Takes no arguments

  --pairs  Compute coverage from paired alignments counting the entire insert as covered.  When using this option only
           reads marked "proper pairs" are used.


Notes:

The input file must be sorted by start position. The samtools package can
be used to sort .bam files. Other files types can be sorted with the "sort"
command (see below).


Example:
   igvtools count -z 5 -w 25 -e 250 alignments.bam  alignments.cov.tdf  hg18

---------------------------------------------------------------------------
Command "sort"
---------------------------------------------------------------------------

Sorts the input file by start position. This command supports the following
file formats:  .bed, .gff, .cn, .igv, .sam, and .bam files


Usage:

  igvtools  sort [options] [inputFile]  [outputFile]

The special string "stdout" can be used as [outputFile], in which case the output will
be written to the standard output stream instead of a file.

Options:

  -t, --tmpDir tmpdir  Specify a temporary working directory.  For large input files
             this directory will be used to store intermediate results of
             the sort. The default is the users temp directory.

  -m, --maxRecords number  The maximum number of records to keep in memory during the
             sort.  The default value is 500000.  Increase this number
             if you receive "too many open files" errors.   Decrease it
             if you experience "out of memory" errors.


---------------------------------------------------------------------------
Command "index"
---------------------------------------------------------------------------

Creates an index for an alignment or the bed feature file formats.  Indexes
required for loading alignment files into IGV, and can significantly
improve performance for large feature files. The input file must be
sorted by start position.  This command does not take an output file
argument, rather the filename is generated by appending ".sai" (for alignments)
or ".idx" (for features) to the input filename. IGV relies on this naming
convention to find the index.

Supported file formats include .bed, .gff, .vcf, .sam., and .bam files


Usage:

  igvtools index [inputFile]


---------------------------------------------------------------------------
Command "formatexp"
---------------------------------------------------------------------------

Format GCT or RES files for display. This should only be used if the file has not previously been log-transformed and has no negative numbers. The module:

1. Takes the log2 of the data.
2. Computes the median and subtracts it from each log2 probe value (i.e., centers on the median).
3. Computer the MAD (mean absolute deviation) using the definition here: http://stat.ethz.ch/R-manual/R-devel/library/stats/html/mad.html
4. Divides each log2 probe value by the MAD.

Supported input file formats are: .gct and .res

Usage:
    
    igvtools formatexp [inputFile] [outputFile]
    
---------------------------------------------------------------------------
Command "gui"
---------------------------------------------------------------------------

Start the igvtools gui

Usage:
    
    igvtools gui
    
---------------------------------------------------------------------------
Command "help"
---------------------------------------------------------------------------

"igvtools help" will display a list of available commands. "igvtools help [command]"
displays help on a particular command.

Example:
    
    igvtools help index

 ---------------------------------------------------------------------------
 Command "version"
 ---------------------------------------------------------------------------

  Prints the igvtools version number.








