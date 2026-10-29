**Track Hubs in IGV**

A track hub is a collection of tracks hosted on a web server and described by a `hub.txt` file, in the format defined
by UCSC (https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html).  IGV can load tracks from a hub, and can also
load a genome a hub defines.  The URL to the `hub.txt` file is the hub URL you give IGV.

**Hubs and genomes**

A hub lists the genomes it provides tracks for.  Each genome is one of:

- A genome already known to IGV, named by its UCSC id (for example `hg38`) or its GenBank/RefSeq accession (for
  example `GCF_000001405.40`).  
- A custom assembly the hub defines itself, with its own sequence (a `twoBitPath`).  Loading one of these creates a
  new genome in IGV, and the hub becomes that genome's annotation hub, supplying its default annotation tracks.

A single hub can list one genome or several, and can mix the two kinds.

**Loading a hub**

Enter the hub URL from either:

- **Tracks > Add Track Hub from URL...**
- **Genomes > Load Genome from URL...**


What happens next depends on the hub.

*Hub for a single genome*

- If the hub's genome is the current genome, the hub's track selection dialog opens.  Check the tracks to load and
  press OK.
- Otherwise IGV switches to the hub's genome, creating it if the hub defines it with a twoBitPath, or loading it from
  the hosted genome list if known to IGV.  Switching genomes unloads all current tracks, so if the hub is loaded
  from the Tracks menu IGV asks for confirmation first.  From the Genomes menu the switch is expected and happens
  without a prompt.  After the switch the hub's track selection dialog opens; the genome's own default annotation
  prompt is not shown here, even for a genome loaded for the first time, so there is one chooser rather than two.
  That prompt waits for the first time the genome is picked from the genome pulldown.

*Hub for several genomes*

- IGV lists the hub's genomes that are not yet in your genome list and offers to add them.  All are checked by
  default; uncheck any you do not want.  Genome ids IGV does not recognize are noted in the dialog.  Cancel adds
  nothing and stops the load.
- Added genomes appear in the genome pulldown.  The first time you load one from the pulldown, IGV prompts for its
  default annotation tracks, as it does for genomes loaded from the host.  That prompt is skipped when restoring a
  session, and if the hub is unreachable the genome loads without annotations and IGV says so.  Annotation tracks
  can be chosen later from the Genomes menu.
- If the current genome is one of the hub's genomes, the hub's track selection dialog opens.  Otherwise IGV offers
  to switch to one of the hub's genomes, or you can cancel and stay where you are.  Switching skips the genome's
  default annotation prompt and opens the hub's track selection dialog instead; the prompt waits for the first load
  from the pulldown.

**The Tracks menu**

Once a hub is known to IGV it appears in the **Tracks** menu, under the **Hubs** heading, whenever the hub defines
tracks for the current genome.  Selecting it opens the hub's track selection dialog.  The info icon beside a hub opens the
hub's description page in your browser, when the hub provides one.

In the track selection dialog, tracks are grouped as the hub organizes them.  Checking a track loads it when you press
OK.  Tracks already loaded are shown checked; unchecking one unloads it.  Cancel leaves everything as it was.

**Select Public Hub**

**Tracks > Select Public Hub...** lists the UCSC public hub registry.  Hubs for the current genome are listed first,
and the list is filtered to the current genome by default; uncheck **Current genome only** to see hubs for other
assemblies.  The filter box narrows the list by name, description or genome.  Hubs already added are greyed out.
Pick one hub and press OK: it is loaded exactly as if its URL had been entered in **Add Track Hub from URL...** (the
hub is added to the Tracks menu and its track selection dialog opens, or IGV offers to switch to the hub's genome).
If you decline the switch you stay on the current genome, but the hub is still added: a hub that defines its own
assembly puts that genome in the genome pulldown, ready to load later; any other hub appears in the Tracks menu once
its genome is loaded.

**Remove Hubs**

**Tracks > Remove Hubs...** lists the hubs you have added, whether selected from the public registry or added by URL,
current genome first.  Check the hubs to drop and press Remove; each is removed from the Tracks menu for every genome
it lists.  Your hubs are saved across sessions in `hubs.txt` in the IGV directory.

A hub that defines a custom assembly is part of that genome rather than managed here, so it is not listed in
Remove Hubs; it is removed with the genome (**Genomes > Remove Genomes...**).  Both menu items are always present;
if there is nothing to list, IGV says so.

**Supported track types**

IGV loads these trackDb types from a hub: bigBed, bigWig, bigGenePred, bigPsl, bigDbSnp, bed, bedGraph, wig,
narrowPeak, broadPeak, gappedPeak, genePred, refGene, gff, gff3, gtf, vcf, vcfTabix, vcfPhasedTrio, bam, tdf, seg,
mut, gistic, interact, rmask and ucscSnp.  Tracks of other types are omitted from the selection dialog.  A hub with no
supported tracks for the current genome reports that rather than opening the dialog.
