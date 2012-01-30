This package contains command line utilities for preprocessing, computing
feature count density (coverage),  sorting, and indexing data files.

---------------------------------------------------------------------------
Starting with shell scripts
---------------------------------------------------------------------------
The utilities are invoked from one of the following scripts:

   igvtools (command line version for linux and  Mac OS 10.x)
   igvtools_gui (gui version for linux and  Mac OS 10.x)
   igvtools_gui.command (alternative double-clickable gui version for Mac OS 10.x)

   igvtools.bat (command line version for windows)
   igvtools_gui.bat (gui version for windows)

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
allocated

   java -Xmx1500m  -jar igvtools.jar [command] [options][arguments]

To start with a gui the command is

   java -Xmx1500m  -jar igvtools.jar -g
   
---------------------------------------------------------------------------
Memory settings
---------------------------------------------------------------------------

The scripts above allocate a fixed amount of memory.  If this amount is not
available on your platform you will get an obscure error along the lines of
"Could not start the Virtual Machine".   If this happens you will need to
edit the scripts to reduce the amount of memory requested,  or use the java
startup option.  The memory is set via a "-Xmx" parameter. For example
-Xmx1500m  requests 1500 MB,  -Xmx1g requests 1 gigabyte.


---------------------------------------------------------------------------
COMMANDS
---------------------------------------------------------------------------

The recognized commands are tile, count, sort, and index.  Note that these
utilities are for working with ascii file formats, including SAM, but
do not work with BAM files.  For manipulating BAM files use samtools.

---------------------------------------------------------------------------
Command "tile"
---------------------------------------------------------------------------

The "tile" command converts a sorted data input file to  a binary tiled
data (.tdf) file. Input file formats supported  are .wig, .cn, .igv,
and .gct.

Usage:

  igvtools tile [options]  [inputFile] [outputFile] [genome]


Required arguments:

  inputFile    The input file (see supported formats below).

  outputFile   Binary output file.  Must end in ".tdf".

  genome       A genome id or filename. See details below. Default is hg18.

Options:


  -z num       Specifies the maximum zoom level to precompute. The default
               value is 7 and is sufficient for most files. To reduce file
               size at the expense of IGV performance this value can be
               reduced.

  -f  list     A comma delimited list specifying window functions to use
               when reducing the data to precomputed tiles.   Possible
               values are  min, max,  mean, median, p2, p10, p90, and p98.
               The "p" values represent percentile, so p2=2nd percentile,
               etc.

  -p file      Specifies a "bed" file to be used to map probe identifiers
               to locations.  This option is useful when preprocessing gct
               files.  The bed file should contain 4 columns:
                  chr start end name
               where name is the probe name in the gct file.


Example:

      igvtoolsh tile -z 5  copyNumberFile.cn copyNumberFile.tdf hg18


Notes:

Data file formats, with the exception of .gct files, must be sorted by
start position.  If neccessary files can be sorted with the "sort" command
described below.  Attempting to preprocess an unsorted file will result
in an  error.

---------------------------------------------------------------------------
Command "count"
---------------------------------------------------------------------------

The "count" command computes average feature density over a specified
window size across the genome. Common usages include computing coverage
for alignment files and counting hits in Chip-seq experiments. Supported
file formats are .sam,  .bam,  .aligned,  .sorted.txt,  and .bed.

Usage:

  igvtools count [options] [inputFile] [outputFile] [genome]

Required rguments:

  inputFile    The input file (see supported formats above).

  outputFile   Binary output file.  Must end in ".tdf".

  genome       A genome id or filename. See details below. Default is hg18.

Options:

  -z num       Specifies the maximum zoom level to precompute.

  -w num       The window size over which coverage is averaged. Defaults
               to 25 bp.

  -e num       The read or feature is extended by the specified distance
               in bp prior to counting. This option is useful for chip-seq
               and rna-seq applications. The value is generally set to the
               average fragment length of the library.

  -f  list     A comma delimited list specifying window functions to use
               when reducing the data to precomputed tiles.   Possible
               values are  min, max,  mean, median, p2, p10, p90, and p98.
               The "p" values represent percentile, so p2=2nd percentile,
               etc.


Notes:

The input file must be sorted by start position. The samtools package can
be used to sort .bam files. Other files types can be sorted with the "sort"
command (see below).


Example:
   igvtools count -z 5 -w 25 -e 250 alignments.bam  alignments.cov.tdf  hg18

---------------------------------------------------------------------------
Genome
---------------------------------------------------------------------------

The genome argument in the tile and count command can be either an id, or
a full path to an IGV .genome file.  The id for IGV supplied genomes are
listed below.  Genome definitions corresponding to these files are in the
"genomes" subdirectory of the igvtools install.  The id is derived by removing
the .extension from the filename.


---------------------------------------------------------------------------
Command "sort"
---------------------------------------------------------------------------

Sorts the input file by start position. This command supports the following
file formats:  .cn, .igv, .sam, .aligned, and .bed.

NOTE: This command will not sort a binary (BAM) file.  Use samtools to sort
and index BAM files.


Usage:

  igvtools  sort [options] [inputFile]  [outputFile]


Options:

  -t tmpdir  Specify a temporary working directory.  For large input files
             this directory will be used to store intermediate results of
             the sort. The default is the users temp directory.

  -m number  The maximum number of records to keep in memory during the
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

Supported file formats are .sam, .aligned, .sorted.txt,  and .bed.


NOTE: This command will not index a binary (BAM) file.  Use samtools to sort
and index BAM files.

Usage:

  igvtools index [inputFile]



 ---------------------------------------------------------------------------
 Command "version"
 ---------------------------------------------------------------------------

  Prints the version to the console.








