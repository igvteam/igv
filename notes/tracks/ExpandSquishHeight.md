**Track height on Expand, Squish and Collapse**

Applies to tracks with rows: feature (annotation) tracks, alignment tracks, variant tracks with genotypes, and
segmented data tracks. Introduced in PR #1901.

***What changed***

Previously, changing the display mode of a track kept the track at the same height. Expanding a track changed
only its content, so the extra rows were hidden behind a scrollbar, and squishing or collapsing a track left
empty space below the content.

Now the track height follows the display mode:

- **Squish** and **Collapse** shrink the track to fit its content. A track never grows on Squish or Collapse,
  so a squished track that still does not fit keeps its height and its scrollbar.
- **Expand** grows the track back to the height it had before it was squished or collapsed. Annotation tracks
  also grow to show up to five rows of features, even if they were never squished. A track never shrinks on
  Expand.
- If an annotation track has no features loaded when you expand it, for example because you are zoomed out
  beyond its visibility window, the track grows when features next load.

These rules apply to the display mode items in a track's popup menu, to the global squish toggle, and to the
`squish`, `expand` and `collapse` batch commands.

***Explicit heights***

A height you set through **Set Track Height...** in the popup menu, through the menu item that sets the height
of all tracks, or with the `setTrackHeight` batch command is kept. The display mode no longer changes the height
of such a track. Setting a height of 0 returns the track to automatic sizing.

Resizing a track by dragging the divider below it does not count as setting an explicit height. A dragged
track still shrinks on Squish and grows on Expand.

***Sessions***

Sessions store track heights as before. Whether a height was set explicitly is not saved, so after a session
is loaded every track follows the display mode rules again. Because Expand only grows and Squish only shrinks
a track, this can only adjust a track toward its content.
