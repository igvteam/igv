
**Background**


Track hubs are described in detail at https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html

All hubs must include a top level "hub.txt" file.  The URL to the hub.txt file is the "hub URL" that you provide to IGV.  
If not using the "useOneFile" option, the hub.txt file points to a genomes.txt file, which in turn points to one or
more trackDb.txt files that define the tracks for each genome listed.  If the hub defines tracks for a single
genome the "useOneFile on" option can be set and contents of the genomes.txt and trackDb.txt files can be included
in the hub.txt file.

***Track data hub***

A track hub defines tracks for a genome assembly that UCSC already hosts, such as hg38, mm39, or dm6.


- You supply the data files (bigBed, bigWig, BAM/CRAM, VCF, HiC, and so on) and the configuration files: `hub.txt`, `genomes.txt`, and one or more `trackDb.txt`.
- In `genomes.txt`, each `genome` line names an existing UCSC database, for example `genome hg38`.
- The reference sequence, base annotations, and tools such as BLAT all come from UCSC's own copy of that assembly.

***Assembly hub***

Assembly hubs are used for assemblies not hosted at UCSC.  They are described at https://genome.ucsc.edu/goldenpath/help/assemblyHubHelp.html
An Assmebly hub is similar to a track hub, but it also provides the following information in the `genomes.txt` file
or `genomes` stanza of the `hub.txt` file:

    - `genome`
    - `twoBitPath`
    - `organism`
    - `scientificName`
    - `description`
    - `defaultPos`
    - `chromSizes`
    - `chromAliasBb`
    - `groups`

**In short**

A track hub layers your data on a genome UCSC already has. An assembly hub brings its own genome sequence, so you can browse an organism or assembly UCSC doesn't carry.

Before building an assembly hub, check UCSC's GenArk collection, which already provides ready-made assembly hubs for thousands of NCBI assemblies. If your genome is there, you can attach a plain track hub to that GenArk assembly instead.