
**Background**


Track hubs are described in detail at https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html

All hubs must include a top level "hub.txt" file.  The URL to the hub.txt file is the "hub URL" that you provide to IGV.  
If not using the "useOneFile" option, the hub.txt file points to a genomes.txt file, which in turn points to one or
more trackDb.txt files that define the tracks for each genome listed.  If the hub defines tracks for a single
genome the "useOneFile on" option can be set and contents of the genomes.txt and trackDb.txt files can be included
in the hub.txt file.


***Genome stanza fields required by IGV***

All genome stanzas must include a `genome` field, which becomes the IGV genome `id`.  If the genome is a custom assembly
it must also include a `twoBitPath` field, which becomes the `twoBitURL` in the genome json. The remaining fields are 
optional, but the genome is degraded without them:

- `chromSizes` -- the 2bit reader does not expose a chromosome list, so without this the chromosome pulldown and 
  whole genome view are disabled.  Files over 10 MB are ignored.  (Note the hub conversion sets 
  `wholeGenomeView = false` regardless.)
- `defaultPos` -- initial locus.  Without it IGV opens on the first known chromosome, if any.
- `scientificName`, `organism`, `description` -- used in that order to build the display name as "name (id)".  With 
  none of them the name is just the id.
- `chromAliasBb` -- chromosome alias lookup.  Without it IGV falls back to its built-in alias defaults.
- `twoBitBptURL` (or `twoBitBptUrl`) -- external 2bit index.  The 2bit's internal index is used otherwise.
- `htmlPath` (-> `infoURL`), `blat`, `taxId` -- informational only.

Path-valued keys are resolved relative to the hub.txt or genomes.txt location by `HubParser`, so relative paths are 
fine. 

**UCSC container types and how IGV groups tracks**

The trackDb spec (https://genome.ucsc.edu/goldenpath/help/trackDb/trackDbHub.html) defines several ways of grouping
tracks.  IGV has no runtime counterpart to any of them -- every hub track becomes an ordinary, independent IGV track.
The container relationships are used only to organize the checkboxes in `TrackSelectionDialog`.  The translation
happens in `TrackDbHub.getGroupedTrackConfigurations`, which turns the flat list of stanzas into a tree of
`TrackConfigContainer` objects (name, label, priority, defaultOpen, tracks, children).  The dialog renders that tree.

*How the tree is built*

- `HubParser` ignores indentation.  Hierarchy comes solely from each stanza's `parent` setting; only the first word of
  the value is kept (`parent myComposite on` -> `myComposite`), so the UCSC `[on/off]` default-selection flag is
  discarded.  `Stanza.getParent()` is resolved after all stanzas are read, so forward references work for inheritance.
- A stanza is a *container* (not a track) if it has any of: `superTrack` without its own `bigDataUrl`, 
  `compositeTrack`, `view`, or `container multiWig`.  Anything else with a `bigDataUrl` and a supported `type` is a 
  track; `cytoBandIdeo` and `assembly` are filtered out.  Container names must be unique (duplicate -> RuntimeException).
- Containers are placed in a parent as they are encountered, so a `parent` container must appear *before* its children
  in the trackDb (groups are the exception, see below).  A container whose parent cannot be resolved is promoted to the
  top level.  A track whose parent cannot be resolved falls into the hub-level container.
- After building, empty containers are dropped and top-level containers are sorted by `priority`.  `priority` strings
  like `3.4` become `3*100 + 4*10` (see `getPriority`); a missing priority sorts last.
- `defaultIsClosed 0` on any container marks it open by default.  `CollapsiblePanel` also opens a panel if it holds
  100 or fewer tracks or any of its tracks are already loaded.

*Setting inheritance*

`Stanza.getProperty` walks up the `parent` chain, matching the UCSC "closest to home" rule, with these IGV
specifics: `track`, `type`, `shortLabel`, `longLabel`, `bigDataUrl`, `parent`, `superTrack`, `priority`,
`compositeTrack`, `view`, `compositeContainer` are never inherited; `visibility`, `priority` and `group` are
*parent-overrides* (the parent's value wins over the child's, mirroring UCSC's rule that composite subtrack visibility
lines are ignored); and `noInherit` on a stanza disables all inheritance.

*Type-by-type translation*

| UCSC construct | Declared by | UCSC meaning | IGV translation |
|---|---|---|---|
| **Hub** (hub.txt `shortLabel`/`longLabel`) | hub stanza | The hub itself; not a track container | A synthetic top-level `TrackConfigContainer` labeled with the hub's `longLabel` (if < 50 chars, else `shortLabel`), priority 0, open.  Holds every track that has no resolvable parent and no group.  Always first in the dialog. |
| **Group** (`groups.txt` stanzas: `name`, `label`, `priority`, `defaultIsClosed`; tracks opt in with `group <name>`) | genome stanza `groups` + track `group` | Browser-wide sections under the image ("Genes and Gene Prediction", etc.).  Hub tracks normally all land in one hub group. | Each `groups.txt` stanza becomes a top-level container (`label`, `priority`, `defaultIsClosed`).  A track or container with no `parent` is assigned to its `group`.  A `group` value not defined in `groups.txt` (or used by a hub with no `groups.txt` at all) creates an implicit container with label = name, priority 1000, open by default.  The raw value is also copied to `TrackConfig.group`. |
| **Supertrack** (`superTrack on [show]`; children `parent <superTrack>`) | parent stanza + child `parent` | A folder, closed by default unless `show`; may contain composites and multiWigs.  One level only. | A container.  Top-level unless it has a `parent` or `group`.  `show` is ignored; open/closed follows `defaultIsClosed` and the panel rules above.  A stanza with `superTrack` *and* its own `bigDataUrl` is treated as a plain track. |
| **Composite** (`compositeTrack on`; subtracks `parent <composite> [on/off]`) | parent stanza + child `parent` | Set of same-type subtracks sharing settings and one configuration page; `[on]` marks a subtrack checked by default. | A container nested under its supertrack/group, or top-level if it has neither.  The `[on/off]` flag is discarded (see *Presentation* below for what actually drives the initial checked state). |
| **Subgroups / dimensions** (`subGroup1..9`, `subGroups`, `dimensions`, `sortOrder`, `filterComposite`) | composite + subtrack | The checkbox matrix on the composite configuration page | Not modeled.  Subtracks are listed flat under the composite.  (The separate `metadata` setting *is* parsed into `TrackConfig.attributes` and shown as columns in the dialog's Search table.) |
| **View** (`subGroup1 view ...` plus a stanza with `view <tag>` and `parent <composite>`; subtracks `parent <viewName>`) | view stanza + child `parent` | A typed sub-division of a composite (e.g. Peaks vs Signals) carrying its own settings | A container nested under the composite.  The `subGroup1 view` declaration is ignored; only the view stanza matters. |
| **Faceted composite** (`compositeTrack faceted`, `metaDataUrl`, `dataTypes`) | composite | Subtracks organized by facets from a TSV file | Treated as an ordinary composite container.  The TSV is not read. |
| **multiWig** (`container multiWig`; subtracks `parent <container>`; `aggregate`) | parent stanza + child `parent` | Several bigWigs rendered as one overlay/stacked track | A container only.  Each subtrack is a separate, individually selectable bigWig; `aggregate` is ignored and no overlay or merged track is created.  `autoScale group` on a subtrack sets `TrackConfig.autoscaleGroup` to the parent's name, so group autoscale does survive. |
| **Predefined groups** (`group map`, `group genes`, ...) | track | Built-in browser groups, "NOT FOR HUBS" | Same handling as any `group` value (see Group row). |

*Presentation in TrackSelectionDialog*

- Each top-level container is one `CollapsiblePanel`: a header reading "label   (N tracks, M selected)" with a
  +/- toggle, followed by the content.
- Inside a panel, each container that directly holds tracks gets a titled-border sub-panel of wrapped checkboxes.
  The title is a breadcrumb built from the ancestor labels, e.g. `My Supertrack - My Composite - Signals`.  Nesting
  depth is conveyed only by the breadcrumb; sub-panels are laid out in a single vertical column, in trackDb order.
- A container with no direct tracks (a supertrack holding only composites) produces no sub-panel of its own; it only
  appears as the first segment of its descendants' titles.
- Checkbox state is `TrackConfig.visible`.  `TrackDbHub` initializes it from the inherited `visibility` setting, but
  `resetSelectionBoxes` overwrites it whenever the dialog is shown: a box is checked if its URL is already loaded,
  or (nothing loaded, initial genome load) if the panel is the `genes` group and the track id is one of `ncbiRefSeq`,
  `augustus`, `ensGene`, or the panel holds a single track.  So in practice `visibility` does not drive the initial
  selection.  The Search button flattens everything into a `TrackChooser` table with a "Group" column equal to the
  top-level panel label.

*What survives into the loaded track*

`TrackConfig.stanzaParent` records the name of the track's top-most ancestor stanza and `TrackConfig.group` records
the `group` setting; both are serialized to the genome/session JSON but are not otherwise used.  Beyond that and the
autoscale group, no container information reaches the `Track` objects.
