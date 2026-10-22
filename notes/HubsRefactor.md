**Track hub refactor**

***Background***

This document describes the refactor of track hub loading to support hubs that reference multiple genomes, and to remove 
the notion of an "assembly hub".  At the end of the refactor there will be a single concept, "Track Hub", which can
reference one or more genomes, each of which can be either a UCSC hosted genome (including GenArk) or a custom
assembly.  The refactor will also remove the `Hub.isAssemblyHub()` method and all logic associated with it.  We will
consider two hub types, (1) a hub referencing a single genome, and (2) a hub referencing multiple genomes.  

***Track hub loading***

A track hub is loaded by URL from the menu, currently from either “File > Load from URL…”, “Genome > Load from URL…”, or
“Track Hubs > Add Track Hub from URL…”. Having 3 items for the same action is possibly confusing and we might
consolidate these, but all 3 items route through “LoadFromURLMenuAction > loadTrackHub” so we will begin our work there.
Leave menu organization for a future task.

NOTE: in this document "assembly" and "genome" are used interchangeably.

A hub can reference a single or multiple assemblies in genomes.txt, or the genome stanza of a useOneFile hub. These
assemblies can be either (1) UCSC hosted genomes, in which case they are referenced by UCSC database id (e.g. hg38) or
genbank accession (e.g. GCA_054883195.1), or (2)  custom assemblies. Custom assemblies require a twoBitPath, and
generally have additional fields such as scientific name. NOTE: A single hub can have a mix of genome types (
see https://bds.mpi-cbg.de/hillerlab/CESAR/trackHub/genomes.txt). 

Currently IGV assumes a hub is either an assembly or a track hub exclusively. This assumption is flawed and we need to 
refactor to account for mixed types.  Consequently we will no longer designate a hub as an assembly hub, the method 
“Hub.isAssemblyHub()”  will go away and all logic associated with it. 

We will now consider 2 cases, (1) a hub referencing a single genome, and (2) a hub referencing multiple genomes.

***Single genome***

The genome might be referenced by ID / Accession or as a custom assembly with a twoBitPath and other fields. IGV already
has code to handle these cases that we can use. All genome stanzas must contain a “genome” field, and will usually have
a trackDb file reference.

If the “genome” property matches the current IGV genome ID we present the track selection dialog from this hub. Also, we 
update the local hub registry to associate the hub with this genome for future sessions
add an entry to the “Track Hubs” menu for this hub. The menu entry should be in current code, the hub registry update
possibly not.

If the “genome” property does not match the current IGV genome we will switch to the hub genome. How this happens
depends on the hub genome type, specifically if it has a “twoBitPath” or not. If it has a “twoBitPath” it is what was
previously called an “assembly” hub and we load it as IGV does currently. If it is lacking “twoBitPath” we assume the
genome property of the hub is either a UCSC genome ID or a genbank/refseq accession known to GenArk. We will refer to
these as genomes “known to IGV”. We already present the complete list of genomes known to IGV in the “Genomes > Load
Genome from Host…” menu. If the genome is in this list we simply load it as if we had selected it there. Additionally,
we update the local hub registry to associate this hub with the current genome.

***Multiple genomes***

If the hub references multiple genomes we will build a list of genomes referenced and present this list to the user to
select genomes they wish to add to the Genomes pulldown. This list will contain any genomes referenced that are not
already in the list. We will build this list from information in the hub itself. For UCSC hosted genomes, including
GenArk, we will probably only have its id (the hubs “genome” property). For genomes with twoBitPath we will usually get
other descriptive fields. Use these fields as we do when generating a name from a hub genome description in
HubGenomLoader. The selection list should present a “select all / none” option, and all should be checked by default.
The dialog would say something like “This hub references the following genomes, add them to the genome selection list?”.
On dismissal IGV would either fetch the genome json for ID/Acession entries or create genome json for “twoBitPath”
entries for the local genome json files, as if they had been loaded from the “Load Genome From Host…” action. We will
not, however, open the track selection dialog for choosing default annotations. Instead, the genome.json will be tagged
so that default annotations are selected the first time it is loaded from the genome selection pulldown. As above in the
single genome case, the hub url itself should be added to the local hub registry to associate this hub with the 
referenced genomws. The genome selection list should be rebuilt to include these new genomes.

NOTE: For genome entries without a twoBitPath we need to verify that they are known to IGV. This is done as described
above in the “Single genome” case. Any genome ids not known to IGV should be reported to the user, possibly in the same
selection dialog as a note.

After updating the genome list we will do one of the following.   (1) If one of the genomes referenced in the hub is the
current genome open the track selection dialog for that hub for the current genome.  (2) Otherwise, we will present the
user with the list of hubs from the genome and ask if they want to switch genomes to one in the list. We must also warn
them that doing so will clear any currently loaded tracks. There should also be an option to cancel and do nothing.

***References***

https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html
https://genome.ucsc.edu/goldenpath/help/assemblyHubHelp.html



***Test Hubs***

https://public.gi.ucsc.edu/hausslerlab/t2t-h9-hub/hub.txt
https://bds.mpi-cbg.de/hillerlab/CESAR/trackHub/hub.txt

**Implementation Plan**

***Current state***

- `LoadFromURLMenuAction.loadTrackHub` parses the hub, then branches on `isAssemblyHub()`.  Assembly hubs go to
  `HubGenomeLoader.loadAssemblyHub`.  Anything else is registered with `HubRegistry.addUserHub` for the current genome,
  with no check that the hub actually references that genome, and no track selection dialog (the call is commented out).
- "Genome > Load from URL" with a hub.txt goes through `HubGenomeLoader.loadGenome(url)`, which always calls
  `loadAssemblyHub` and throws on a track hub.
- `Hub.getGenomeConfigs()` builds a `GenomeConfig` for every genome stanza, including id-only stanzas with no sequence.
  A mixed hub therefore produces invalid genome configs.
- `GenomeSelectionDialog.loadButtonActionPerformed` holds the "load from host" logic (fetch config, prompt for
  annotations, `GenomeDownloadUtils.downloadGenome`, load).  It is private and reads the dialog's own radio buttons.
- "Known to IGV" already exists: `HostedGenomes.getGenomeListItem(id)` checks the IGV list first and the GenArk list
  second.  `GenomeListManager.getGenomeTableRecord(id)` covers genomes already in the pulldown.  Both lists are already
  loaded by the time the genome selection dialog has been used; in any event the fetch is required.
- `HubRegistry` persists user-selected and user-added hubs in `~/igv/hubs.txt`, tab-delimited:
  `hubUrl, shortLabel, longLabel, dbList, descriptionUrl`.  `dbList` is the comma-separated list of genome ids from
  the hub's genome stanzas.  `IGVMenuBar.updateTracksMenu` shows a registry hub for the current genome when
  `genome.getUCSCId()` matches an entry in `dbList` exactly.  Registry entries are independent of the genome json, so
  they survive a genome refresh from the host.
- Nothing in `GenomeConfig` or `GenomeManager` supports "select default annotations on first load".

***Genome ID aliasing -- review***

Hub genome stanzas name genomes by UCSC database name (hg38) or GenArk accession (GCA_/GCF_).  A hub naming GRCh38 is
not expected and is probably not legal.  The aliasing problem is on the genome side: the current IGV genome may carry a
non-UCSC id (GRCh38 from a user json or session, hg38_1kg, b37, 1kg_v37) and must still be recognized as the hub's
hg38 / hg19.  Current state:

- `Genome.getUCSCId()` returns the json `ucscID` if set, else `Genome.ucsdIDMap` (1kg_ref -> hg18, 1kg_v37 -> hg19,
  b37 -> hg19), else the id.  The hosted `hg38_1kg.json` sets `"ucscID": "hg38"`; no other IGV-hosted genome sets one,
  and all other hosted ids are plain UCSC database names.  So matching the hub's `genome` against `getUCSCId()`
  (never `getId()`) already handles the 1kg case without forcing a genome switch.  `ucsdIDMap` lacks the GRC names.
- `Hub.idMappings` maps UCSC database name -> RefSeq/GenBank accession and is used in `Hub.getTrackDbHub` so an IGV
  genome "hg38" finds a hub stanza named "GCF_000001405.40".  It is incomplete and has suspect entries:
    - missing hosted genomes: hg19, hs1, mm9, danRer10/11, ce10, dm3, bosTau8, canFam3/5, panTro4/5, gorGor4,
      panPan2, strPur2
    - `tetNig2` -> GCF_000002275.3 duplicates the `taeGut1` accession; `xenTro9` -> GCF_000002035.6 duplicates
      `xenTro10`; `sacCer2` -> GCF_000146045.2 duplicates `sacCer3`.  These look like copy errors.
- `HubRegistry` matches `dbList` entries against `getUCSCId()` exactly.  This is correct once `getUCSCId()`
  canonicalizes properly; the registry itself needs no aliasing.
- `EncodeTrackChooserFactory.getEncodeGenomeID` and `Mutation` carry their own hard-coded equivalences.  Out of scope.

Plan:

- Extend `Genome.ucsdIDMap` with GRC-style names so `getUCSCId()` canonicalizes them: GRCh38 -> hg38, GRCh37 -> hg19,
  GRCm38 -> mm10, GRCm39 -> mm39, T2T-CHM13v2.0 -> hs1.  This is the only genome-side change.
- Fix and complete `Hub.idMappings` (database name -> accession) for the UCSC genomes IGV hosts: correct the three
  suspect entries and add the missing genomes.  Each accession to be verified against UCSC/GenArk before committing.
- Matching rule for "does this hub reference the current genome": the hub stanza's `genome` equals
  `currentGenome.getUCSCId()`, or equals `idMappings.get(currentGenome.getUCSCId())`.  No new utility class.


***Steps***

1. **Per-genome classification in `Hub`.**  Add accessors for the genome ids referenced, whether a given id has a
   `twoBitPath`, and `getGenomeConfig(id)` that builds a config only for twoBit stanzas.  Delete `isAssemblyHub` and
   `HubParser.loadAssemblyHub`.  Use the matching rule from the aliasing section.

2. **Extract the host-load logic.**  Move the body of `GenomeSelectionDialog.loadButtonActionPerformed` into a static
   method, roughly `loadHostedGenome(GenomeListItem, downloadSequence, downloadAnnotations, selectAnnotations)`, that
   produces the local genome json and optionally loads it.  The dialog calls it unchanged.

3. **Recording the hub.**  Use `HubRegistry.addUserHub`, not the genome json.  The registry entry's `dbList` carries all
   of the hub's genome ids, so the hub appears in the "Track Hubs" menu for every genome it references, and it survives
   a genome refresh from the host.  `addUserHub` should skip a URL already present (it currently appends blindly).
   The genome json `hubs` property is written only when IGV creates a new json for a custom (twoBit) assembly from the
   hub, as it does today: the hub URL is `hubs[0]`, making it the genome's annotation hub.  This is creating a json,
   not modifying one; an existing genome's json is never modified.  Genomes loaded from a remote json URL or a fasta
   need no special case: their URL is their id, and the registry entry applies to them like any other genome.

4. **First-load annotation flag.**  Add a boolean to `GenomeConfig`, serialized in the json, meaning "prompt for default
   annotations on first load".  `GenomeManager.loadGenome` honors it after the loader returns: run
   `selectAnnotationTracks` (which uses `hubs[0]`, the genome's own annotation hub), set tracks, clear the flag, save
   the json, then continue to `setCurrentGenome`.  Local-only field; igv.js ignores unknown properties.

5. **New orchestrator**, replacing the static half of `HubGenomeLoader`.  One entry point taking a parsed `Hub`,
   dispatching on genome count.
   - Single genome, equivalent to current: show the hub's track selection dialog via
     `SelectHubTracksAction.selectAndLoadTracks`, register with `HubRegistry.addUserHub`, update menus.
   - Single genome, twoBit: existing assembly path (select annotations, save a new json with the hub URL as
     `hubs[0]`, load).
   - Single genome, id-only: look it up with `GenomeListManager.getGenomeTableRecord` then
     `HostedGenomes.getGenomeListItem`; load through step 2; register the hub.  Unknown id: message to the user.
   - Multiple genomes: build the candidate list, excluding genomes already in the pulldown and flagging unknown ids.
     Show the new checkbox dialog.  For each accepted entry, create (twoBit, with the hub URL as `hubs[0]`) or fetch
     (id-only) the json, set the first-load flag, add a `GenomeListItem`.  Register the hub once.  Then either open track selection for the
     current genome if it is equivalent to one of the hub's genomes, or offer the switch-genome prompt with the
     clear-tracks warning and a cancel.

6. **Two new dialogs** in `org.igv.ucsc.hub`: the multi-genome add list with select all/none and an unknown-ids note,
   and the switch-genome chooser.  Verified offscreen by painting to a BufferedImage.

7. **Rewire entry points.**  `loadTrackHub`, the hub branch of `loadGenomeFromUrl`, and `loadUrls` all call the
   orchestrator.  `GenomeSelectionDialog`'s GenArk branch uses `getGenomeConfig(id)` instead of
   `getGenomeConfigs().get(0)`.  Keep the `HubGenomeLoader` instance loader, since `GenomeLoader.getLoader` and
   sessions still route hub URLs through it.

8. **Cleanup and docs.**  Remove `loadAssemblyHub`, the SwingWorker `loadGenome(url)`, and the private duplicate
   `isHubURL`.  Update `HubTest` to use `loadHub`.  Add a mixed-hub fixture under test data with one test for the
   step 1 classification and one for alias equivalence.  Update `notes/HUBS.md` to drop the assembly hub section.

***Resolved: first-load annotation prompt***

The purpose of creating genome json files from the hub's genomes.txt is to populate the genome selection pulldown.
Those json files are not complete until the user has selected initial annotations, which happens lazily on first
load.  Annotation prompting must then work exactly as "Load Genome from Host": the prompt presents the genome's own
annotation hub (`hubs[0]`).  For a custom (twoBit) assembly that is the hub itself; for an id-only genome it is the
hosted definition's hub.  Tracks from the loaded hub are reached through the "Track Hubs" menu.


***Implementation notes (Oct 21 2026)***

Implemented as planned, with these specifics:

- The first-load annotation flag (`selectAnnotationsOnLoad`) is honored in `JsonGenomeLoader.loadGenome`, before the
  `Genome` is constructed, rather than in `GenomeManager.loadGenome` after it: the Genome derives its annotation
  tracks from the config at construction, so the prompt must precede it.
- `Hub.idMappings` was verified against GenArk by fetching each candidate hub.txt.  No GenArk hub exists for hg19,
  hg18, mm9, dm3, danRer10, ce10, canFam5, panTro4/5, gorGor4, taeGut1, xenTro9, so they are not mapped.  hs1 maps
  to GCA_009914755.4 (the GCF accession has no hub).  The pre-existing hg38 entry is kept.  tetNig2 and xenTro10 were
  corrected; sacCer2 was dropped.
- Single genome, id-only, already in the genome list: loaded with `GenomeManager.loadGenomeById`, as selecting it from
  the pulldown would.  Only genomes not yet in the list go through `GenomeSelectionDialog.createLocalGenome`.
- `IGVMenuBar.updateTracksMenu` omits a registry hub that is also one of the genome's own hubs, so an assembly created
  from a hub lists that hub once.
- `HubParser` keeps only the first whitespace-delimited word of most stanza values, so a multi-word `scientificName`,
  `organism` or `description` is truncated in the generated genome name ("Python" for "Python bivittatus").  Pre-existing,
  not changed.
- After switching to a hosted (id-only) genome -- single genome hub, or the switch offered after a multi-genome hub --
  the hub's own track selection dialog is shown, following the first-load default annotation prompt when there is
  one, so the hub's tracks are not overlooked.  Not done for an assembly created from the hub, whose default
  annotations are the hub's tracks.
