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

https://hgdownload.soe.ucsc.edu/hubs/GCA/054/883/195/GCA_054883195.1/hub.txt

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


***Follow-up changes (Oct 22 2026)***

Changes after review and manual testing, all on branch `trackhubs` (PR #1887):

- **Tracks menu.**  The "Hubs" menu was renamed `tracksMenu`; hubs are only part of it.  Registry hubs for the current
  genome are listed under a bold "Hubs" heading, each with an info link-out drawn from the hub's `descriptionUrl`
  (`HubMenuItem`, same icon as the hub selection dialog).  A separator always precedes "Manage Track Hubs ...", which
  is shown whenever any hub is known.  On OK a dialog summarizes the hubs added and removed; nothing is shown if the
  selection is unchanged.
- **Manage Track Hubs lists every hub**, for all genomes, sorted so hubs for the current genome come first.  Filtering
  to the current genome and then saving the result replaced the whole of `hubs.txt`, dropping hubs selected for other
  genomes.  Only registry hubs (UCSC public list + user added) are managed; a genome's own hubs from its genome json
  are not.
- **Assembly hubs are never registry hubs.**  A genome created from a hub's `twoBitPath` records the hub in its own
  json (`hubs[0]`), so `TrackHubLoader` does not also add it to `HubRegistry`; otherwise it is counted twice.  A
  multi-genome hub is registered only if it references at least one genome it does not define, so a mixed hub such as
  the T2T-H9 hub still reaches the Tracks menu for hg19.  This amends step 3 and step 5 above.
- **Accession aliasing in the registry and the multi-genome loader.**  `HubRegistry` lookups for a genome also match
  its GenArk accession (`Hub.getGenomeAlias`), so a hub whose genomes.txt names GCF_000001405.40 is listed under hg38.
  `TrackHubLoader.loadMultiGenomeHub` matches hub genome ids against the genome list the same way, via the reverse
  lookup `Hub.getUcscIdForAlias`, so an existing hg38 is not offered again as a GenArk genome.
- **Genome switch confirmation.**  A single-genome hub referencing another genome switches silently from the Genomes
  menu but asks first from the Tracks menu (`loadHub(url, confirmGenomeSwitch)`), since the switch unloads all tracks.
- **First-load annotation prompt rules** (`JsonGenomeLoader.selectAnnotations`).  The `selectAnnotationsOnLoad` flag is
  cleared only once the prompt has been shown.  If it cannot be shown (batch, headless, tests, no IGV instance) the flag
  and json are left alone for a later interactive load.  If the annotation hub is unreachable the genome still loads,
  the user is told with the error, and the flag is cleared so they are not asked on every load; annotations can be
  added later from the Genomes menu.  No prompt during session restore, JSON or XML: `GenomeManager` passes a
  `promptForAnnotations` flag to the loader, false from `loadGenomeById(.., loadAnnotationTracks=false)` (JSON
  sessions) and from `loadSessionGenome` (XML sessions, which still load default annotations).
- **Multi-genome dialog flow.**  Cancel in "Add Genomes from Hub" aborts the load entirely.  The dialog is constructed
  and shown on the event thread, and the switch-genome prompt is marshaled there too.
- **Track selection dialog.**  All checkboxes are enabled; unchecking a loaded track unloads it on OK, as the ENCODE
  chooser does.  The diff is a single pass over the hub's configs with hash lookups against loaded track paths.
- **Formats.**  `bigPsl` added to `TrackDbHub.supportedTypes` (the generic BigBed codec handles bed12+13).  The T2T-H9
  hub's hg19 trackDb, reached through an `include`, was being rejected for this reason, not the include.


***Copilot review fixes (Oct 23 2026)***

Six findings from the Copilot review of PR #1887, all confirmed and fixed in one commit:

- XML session restore went through the one-argument `GenomeManager.loadGenome`, so the first-load annotation prompt
  still fired.  Added `loadSessionGenome` (annotations loaded, no prompt) for `XMLSessionReader`.
- The single-genome path's `findKnownGenome` matched the genome list by raw id, so a hub naming hg38 by accession
  loaded a duplicate GenArk genome.  It now uses the alias-aware `findInGenomeList` first.
- `SelectHubTracksAction.selectAndLoadTracks` built and showed the track selection dialog on the calling thread, which
  is a background thread from both the hub loader and the menu action's SwingWorker.  Dialog creation and display are
  now marshaled to the event thread.
- Failed genome hubs were listed inside the block guarded by a non-empty loaded-hub list, so a genome whose only hub
  failed showed nothing.  The block now opens when either list is non-empty.
- Two trackDb entries can share a URL; unchecking one removed the track even if the other stayed checked.  Checked
  URLs are now subtracted from the removal set.
- `TrackSelectionDialog` cached dialogs by hub only.  A mixed hub opened for one genome and then another reused the
  first genome's dialog while the action read the second genome's configs.  The cache is now keyed by hub url and
  genome id.

***Hub package review fixes (Oct 23 2026)***

A second, medium-effort `/code-review` scoped to `org.igv.ucsc.hub` after the Copilot round.  One finding was rejected:
cancelling the default-annotation dialog does not cancel a genome switch; the genome is saved with no tracks and
loaded, since a genome without annotations is unusual but allowed.  The rest, fixed in one commit:

- `Hub.getGenomeConfig` returns the config cached inside the (session-cached) Hub, and three sites edited it in place
  (`hubs`, `tracks`, `selectAnnotationsOnLoad`), so one flow's edits leaked into the next -- e.g. a GenArk genome added
  from a multi-genome hub left the flag set, and loading that GenArk hub later prompted twice.  All three sites
  (`TrackHubLoader` x2, `GenomeSelectionDialog.createLocalGenome`) now take `copy()` first.
- `HostedGenomes.getGenomeListItem` keys IGV's own hosted list by UCSC name, so a hub naming hg38 by accession fell
  through to the GenArk record.  `TrackHubLoader.findHostedGenome` tries the UCSC name (`Hub.getUcscIdForAlias`) first.
  Hosted candidates now carry the resolved id so the genome list entry matches the json.
- `HubMenuItem`: `setToolTipText(null)` (a hub with no longLabel) unregistered the item from the ToolTipManager and
  silenced the info icon's url tooltip.  The item re-registers itself when it has an info url.
- "Add Genomes from Hub" opened with an empty list when every referenced genome was already in the list and only
  unknown ids remained, and Cancel then aborted the load.  The unknown ids are now reported in a message and the load
  continues.
- `Hub.getTrackDbHub` cached the parsed trackDb under the caller's id, so lookups by UCSC name and by accession each
  fetched and parsed it.  It is cached under the stanza's own `genome` value, which both lookups reach.
- The hosted-genome branch of `loadSingleGenomeHub` registered the hub before loading the genome; a failed load left it
  in `hubs.txt`.  Registration now follows the load.

***GRC alias gaps (Oct 23 2026)***

The GRC <-> UCSC aliases (GRCh38 -> hg38 etc.) are the ones that matter in practice; accession aliases are rarely hit.
The current-genome path already canonicalized through `Genome.getUCSCId` (covered by `HubTest`).  Two edges did not:

- `GenomeConfig.getUcscID` returned the raw id, so a json genome with id GRCh38 and a `hubs` entry could not find its
  hg38 trackDb in the annotation prompt.  The map is now behind a static `Genome.toUCSCId`, used by both classes.
- `TrackHubLoader.findInGenomeList` matched list entries by raw id (plus accession), so a list entry named GRCh38 or
  b37 was not recognized as a hub's hg38 / hg19 and the hosted genome was offered again.  It now compares canonical
  UCSC names across the list.

Hubs themselves are still expected to name genomes by UCSC database name or accession, never GRC names; nothing on the
hub side maps those.

***Second Copilot review (Oct 23 2026)***

Eleven items; two were already fixed (dialog cache per genome, XML session prompt) and only their threads were open.
Fixed in one commit:

- The hosted-genome branch registers the hub after the load, but the load had already rebuilt the Tracks menu, so the
  hub was missing until the next rebuild.  Menus are refreshed after registration.
- `TrackSelectionDialog` starts canceled and OK clears it; a reused dialog kept the previous OK, so closing with the
  window X after an earlier OK unloaded unchecked tracks.  Reuse resets it to canceled.
- The JSON session path tied the annotation prompt to the session's own load-default-tracks option (a regression from
  the first round).  The session overload of `GenomeManager.loadGenomeById` never prompts; the pulldown path still does.
- `GenomeManager.selectAnnotationTracks` and `HubGenomeLoader.selectAnnotationTracks` built the track dialog off the
  event thread with only `setVisible` marshaled.  Both now build and show it on the event thread.
- A genome created from GenArk keeps the accession as its UCSC id, and `Hub.getGenomeStanza` / `HubRegistry` only
  mapped UCSC name -> accession, so an hg38 hub did not match it and the user got a needless switch prompt.  Both now
  fall back to the reverse mapping.  The trackDb cache is keyed by the stanza's name on every path.
- `Hub.getTrackDbHub` swallowed a trackDb fetch failure, so with hub.txt cached but trackDb unreachable the annotation
  chooser opened empty and OK saved no tracks and cleared the flag.  New `Hub.ensureTrackDbLoaded` propagates the
  IOException; the annotation prompt calls it first, so `JsonGenomeLoader` reports an unavailable hub.  The menu
  lookups still swallow it as before.
- "Is this the genome's own hub" was decided by `hasTwoBit`, so a hub defining an assembly with the current genome's id
  (but not the hub the genome was created from) was neither registered nor in the json and vanished after the
  immediate chooser.  Jim's call: test by hub URL against `genome.getTrackHubs()`.  Single-genome: register unless the
  hub is one of the current genome's own hubs.  Multi-genome: register when the hub references any genome that already
  existed (not created by this load), unless it is the current genome's own hub; an all-assembly hub whose genomes
  were all just created stays out.  Consequence: reloading such a hub once its genomes exist registers it and it
  appears in Manage Track Hubs; the Tracks menu still lists it once (registry hubs that are also the genome's own are
  dropped there).

***Manage Track Hubs replaced (Oct 28 2026)***

The multi-select "Manage Track Hubs ..." dialog (check to add, uncheck to remove, over every hub for every genome)
proved confusing and unworkable.  Replaced by two explicit Tracks menu items:

- **Select Public Track Hub ...** -- `HubSelectionDialog` repurposed as a single-select chooser over the UCSC public list
  only (`HubRegistry.getUcscPublicHubDescriptors`).  `HubTable` gained a single-select mode: the checkbox column is
  hidden and rows whose url is in a disabled set (hubs already in the registry) are greyed out and unselectable;
  the descriptors' "selected" flag stays the checkbox state and is not reused.  A "Current genome
  only" checkbox, on by default, filters the list; off, a hub for another assembly can be chosen.  Jim's call: the
  chosen hub is loaded exactly as "Add Track Hub from URL" loads it (`LoadFromURLMenuAction.loadTrackHub` ->
  `TrackHubLoader.loadHub`), so registration, menu refresh, the track selection dialog and the genome-switch
  confirmation are all the same code path.
- **Remove Hubs ...** -- new `RemoveHubsDialog`, modeled on Remove Genomes: a checkbox `HubTable` over
  `HubRegistry.getSelectedHubs()` (UCSC-selected and URL-added alike), "Remove" writes the remainder with
  `setSelectedHubs`.  Shown only when the registry is non-empty.

The added/removed summary message is gone; each action gives direct feedback.

Follow-up: both items are always shown.  "Remove Hubs" was first shown only when the registry was non-empty, and an
assembly hub selected from the public list goes to the genome's json, not the registry, so after adding the sunflower
sea star hub the item was simply absent.  Jim's call: disappearing menu items are confusing; an empty list is reported
with a message instead (and the message points assembly hubs at Remove Genomes).  Same for "Select Public Track Hub" when
the UCSC list could not be loaded.

Follow-up: declining the genome-switch prompt (single-genome hub for another genome, from the Tracks menu or Select
Public Track Hub) used to discard the hub.  Jim's call: the hub is still added and the user simply stays on the current
genome.  A hub that does not define its genome is registered; it appears in the Tracks menu once its genome is
loaded, and is greyed out in Select Public Track Hub.  An assembly hub is instead added by putting its genome in the genome
list with the first-load annotation flag, as the multi-genome path does (Copilot review of #1892: registering it
greyed it out in the chooser with no genome to load, a dead end).

Follow-up (Oct 29 2026): switching to a hosted genome from a hub load no longer shows the genome's first-load
default annotation prompt.  Jim's call: a hosted genome loaded for the first time this way showed two choosers in a
row, the default annotations and then the hub's tracks, which was confusing.  The hub load passes
`promptForAnnotations=false` through the three-argument `GenomeManager.loadGenomeById` that session restore already
used, and `createLocalGenome(..., selectAnnotations=false)`, so the genome loads with whatever `tracks` its json carries (IGV's
hosted json has defaults; a GenArk genome has none) and `selectAnnotationsOnLoad` stays set for the first load from
the pulldown.  Covers the single-genome hosted path and the switch offered after a multi-genome hub.

Same day, extended to assemblies the hub defines (Jim's call): the single-genome assembly path no longer runs
`HubGenomeLoader.selectAnnotationTracks`; it saves the json with `selectAnnotationsOnLoad` set, as the multi-genome
`addGenome` already did, loads without the prompt, and opens the hub's track selection dialog like the hosted path.
The multi-genome switch lost its twoBit branch, which had suppressed the track dialog for assemblies -- backwards,
since the annotation prompt was the dialog being avoided.  The hub's track dialog is now the one chooser for every
hub load.  After `/simplify`: the single-genome assembly path shares `addAssemblyGenome` with the declined-switch
path (genome list entry with the first-load flag, created only if the genome is not in the list yet, so a repeat
load of the hub no longer rewrites a json whose annotations were already chosen) and loads by id.

Bug (Oct 29 2026): the Weng lab ENCODE hub (https://storage.googleapis.com/gcp.wenglab.org/hubs/dna20/hub.txt) opened
an empty track selection dialog for hg38 with no "filtered out" message.  Every one of its 21,699 data tracks carries
a `view` line (even ones directly under a composite, e.g. `view bigWig`), and `TrackDbHub` took any stanza with
`view` as a view container, so each data track became an empty container and was trimmed; the supported-track count,
which ignores containers, still said 21,699, so the dialog opened, and nothing was counted as unsupported.  Not new
to the refactor.  Fix: a stanza with its own `bigDataUrl` is a data track whatever else it declares (the superTrack
test already said so; now the whole container predicate does).  Fixture `test/data/hubs/view_on_datatrack_trackDb.txt`.
