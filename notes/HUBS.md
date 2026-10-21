
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
