**Track heights: current state and deferred preference work**

***Background***

Written 2026-10-29 while fixing hub bigwig tracks that loaded at 12 px (branch `min-track-height`).  That change
collapsed IGV's three competing "default/minimum" heights into one `IGV.track.height` preference (default 20),
removed `AbstractTrack.DEFAULT_HEIGHT`, the broken `Constants.INITIAL_TRACK_HEIGHT` ("15" used as a key), and the
renderer-level minimum (`DataRenderer` 20, `HeatmapRenderer` 1, a vestige of one-track-per-heatmap-row), and
clamped hub `maxHeightPixels` values to the preference in `TrackDbHub.getTrackConfig`.

The follow-on idea, deferred for now, is to give every track type with a hardcoded default or minimum height its
own preference, as `InteractionTrack` (`IGV.interact.track.height`), `AlignmentTrack` (`ALIGNMENT_TRACK_HEIGHT`),
`GWASTrack` (`GWAS_TRACK_HEIGHT`) and `DataTrack` (`IGV.chart.track.height`) already have.  This note records the
survey so that work can start without redoing it.

***How a track gets its height***

- `AbstractTrack.height` has no initializer.  0 means "auto-size": `getHeight()` returns `getContentHeight()`.
  Nothing clamps `setHeight`.  Sessions (XML and JSON) assign the field directly, bypassing `setHeight`.
- `TrackPanel.getPreferredSize` uses `max(track.getHeight(), track.getContentHeight())`.
- `Track.getMinimumHeight()` is enforced only by the divider drag (`TrackPanelDivider`) and `minimizeHeight()`.
  `AlignmentTrack` overrides the minimum to 50.
- Row-based tracks (feature, alignment, seg, variant) have no initial height at all.  They are auto-sized from
  the row count times the expanded/squished row height.  Their "default height" is really a default row height.
- Hub heights: `TrackDbHub.getTrackConfig` parses `maxHeightPixels max:def:min` into `TrackConfig.maxHeight /
  height / minHeight`, now clamped to `IGV.track.height`.  `TrackProperties(TrackConfig)` drops `maxHeight`.
  Genome json `tracks` entries go through `TrackConfig.fromJSON` and are not clamped.  Track lines in data file
  headers go through `ParsingUtils.parseTrackLine` and are not clamped (deliberate: scope was hubs only).
- `TrackDbHub` and `ParsingUtils` read a 2-value `maxHeightPixels` differently (`max:def` vs legacy `def:min`).
  `TrackDbHub` has no `NumberFormatException` handling.  Pre-existing, unchanged.

***Existing height preferences***

| Key | Default | Used by |
|---|---|---|
| `IGV.track.height` | 20 | minimum height for all tracks; default expanded row height fallback; hub clamp; Set Track Height dialog seed |
| `IGV.chart.track.height` | 40 | `DataTrack` constructor |
| `IGV.interact.track.height` | 250 | `InteractionTrack` constructor |
| `ALIGNMENT_TRACK_HEIGHT` | 300 | `AlignmentTrack` constructor |
| `GWAS_TRACK_HEIGHT` | 200 | `GWASTrack` constructor (hidden preference) |

***Hardcoded heights that could become preferences***

Whole-track starting heights (tracks that are not row-based and call `setHeight` explicitly):

| Literal | Where |
|---|---|
| 50 | `GisticTrack.DEFAULT_HEIGHT`, constructor |
| 50 | `CNFreqTrack` constructor |
| 300 | `BasePairTrack` constructor |
| 30 | `CytobandTrack` constructor |
| 60 | splice junction track, set in `AlignmentTrack` constructor |
| 200 | `TrackLoader.ALIGNMENT_FEATURE_HEIGHT`, psl / bigPsl / bigChain tracks when the locator has no display mode |
| 300 | `VariantTrack.DEFAULT_MAX_HEIGHT`, cap on the auto-sized height |
| 15 | `TrackLoader`, Repeat Masker feature files |
| 40 | `TrackLoader`, allele frequency files (with `PointsRenderer`) |
| 45 | `EmblFeatureTableParser`, EMBL feature tracks |

Expanded / squished row heights (auto-sized row tracks):

| Expanded | Squished | Where |
|---|---|---|
| 35 | 12 | `FeatureTrack.DEFAULT_EXPANDED_HEIGHT` / `DEFAULT_SQUISHED_HEIGHT` |
| 14 | 2 | `AlignmentTrack.DEFAULT_EXPANDED_ROW_HEIGHT` / `DEFAULT_SQUISHED_ROW_HEIGHT` |
| 15 | 2 | `SegTrack.EXPANDED_SAMPLE_HEIGHT` / `SQUISHED_SAMPLE_HEIGHT` |
| 15 | 4 | `VariantTrack.DEFAULT_EXPANDED_GENOTYPE_HEIGHT` / `DEFAULT_SQUISHED_GENOTYPE_HEIGHT` |

Minimum heights:

- `AlignmentTrack` stores 50 twice: the `minHeight` field (used in `getContentHeight`) and the
  `minimumHeight = 50` assignment in the constructor (used by `getMinimumHeight`).  One preference should replace both.

Layout constants that are probably not candidates:

- `VariantTrack.VARIANT_BAND_HEIGHT` 25: the fixed variant-site band.  Render and mouse lookup both read it, so
  whatever happens it must stay a single source.
- `SequenceTrack.SEQUENCE_HEIGHT` 17, `IGVFeatureRenderer.BLOCK_HEIGHT` 14 and friends: drawing geometry.

***Open questions for when this is picked up***

- Which families to expose: whole-track defaults, row heights, the alignment minimum, or some subset.  Eight row
  height preferences is a lot of UI for the Tracks tab.
- Key naming.  Existing keys mix styles (`IGV.chart.track.height`, `ALIGNMENT_TRACK_HEIGHT`, `GWAS_TRACK_HEIGHT`).
  Pick one for the new keys.
- Whether `GWAS_TRACK_HEIGHT` should become visible in the Tracks tab alongside the others.
- Whether the `TrackDbHub` / `ParsingUtils` `maxHeightPixels` parsing differences should be reconciled, and whether
  track-line heights should be clamped like hub heights.
