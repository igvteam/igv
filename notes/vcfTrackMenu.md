
## Changes to the VCF track menu

1. Change the label “Color by” to "Color variant sites by”

2. If the user selects an entry under "INFO field", do not pop up a color dialog. Just use the default color for that field.

3. If the user selects an entry under "INFO field", the checkmark displayed next to the "INFO field" menu item should look the same as the checkmark displayed when "Allele frequency" or "Allele fraction" is selected. It should work with the different themes, included dark mode.

4. Replace the menu items "Allele Frequency Colors..." and "Edit INFO Colors..." with just one menu item: "Edit Color..." which combines the functionality of both items that it replaces, by bringing up the color dialog for the currently selected option, which will be one of the following: Allele frequency with color scale, Allele frequency with bar height, Allele fraction with color scale, Allele fraction with bar height, Currently selected INFO field.  Bar height vs color scale for allele frequency and fraction is determined by the "Show allele frequency / fraction as" option in the Variants preferences tab.

5. Rename the "Edit Color..." dialog title to "<selected option> Colors" where <selected option> is the currently selected option, which will be one of the following: Allele Frequency/Fraction or  INFO Field <info field> , where <info field> is the currently selected INFO field.  The dialog title should change when the user selects a different option in the menu.  Also Rename the "Edit Color..." menu item to "<selected option> Colors..." where <selected option> is the currently selected option, which will be one of the following: Allele Frequency/Fraction or  INFO Field <info field> , where <info field> is the currently selected INFO field.  The menu item should change when the user selects a different option in the menu.

6. In the allele freq/frac + bars color dialog, make the following changes: (1) No value or 0% should always have the same color as Reference. No separate color specification is needed for that case. Currently, it uses track color which isn’t editable. Remove that bar from the dialog. (2) Replace the text at the top of the dialog with the following, but add a line break after every sentence.
"The lower part of each bar is colored in the variant color and its height is in proportion to the allele frequency or fraction. Allele frequency = AF; or GMAF if a variant has no AF. Allele fraction = AC divided by AN. A variant with no value is drawn in the reference color. A gVCF reference block is drawn in gray. Filtered variants are drawn faded."

7. In the allele freq/frac + color scale dialog, make the following change: Replace the text at the top of the dialog with the following, but add a line break after the first sentence. "Variants are colored by allele frequency or fraction, blended on a log scale between stops. Variants with no value take the color of the lowest stop."



## Variant preferences

The View > Preferences > Variants tab currently has two boxed sections labeled: *Variant colors* and *Variant INFO colors*.

This needs to change as follows: 

There should be two boxed sections labeled: 
    
1. Genotype colors 
2. Variant site colors

The *Genotype colors* box should contain the first 4 items in the current *Variant colors* box, i.e. Homozygous reference, Homozygous variant, Heterozygous variant, No call and their corresponding color swatches.  Clicking on the color swatches should behave the same way it does now.

The *Variant site colors* box should be divided by a horizontal line into 2 sections.

Above the line should contain all the remaining items from the current *Variant Colors* box that are not in the new *Genotype colors* box. 
Also change the label for the top item to "Allele frequency / fraction shown as bar height - Reference"
and change the label for second item to "Allele frequency / fraction shown as bar height - Variant".  
All items should behave the same way they do now.

Below the line should contain everything from the current *Variant INFO colors* box, with the following change: Remove the first line of text, and change the word "attribute" in the second line of text to "VCF INFO attribute".  
The table and buttons should behave the same way they do now.


## Implementation plan

Everything above the "IGNORE FOR NOW" section. Menu items 1-5 are in `src/main/java/org/igv/variant/VariantTrackMenuHelper.java`; items 5-7 also touch the three dialogs and `VariantRenderer.java`; the preferences changes are in `preferences.tab`, `PreferencesEditor.java` and `VariantColorSchemePanel.java`. Before starting, run `./gradlew test --tests 'org.igv.variant.*' --tests 'org.igv.prefs.*'` for a baseline, and again afterwards.

### Menu items

**1. Label.** In `getMenuItems` (line 70), change `<b>Color By` to `<b>Color variant sites by`. "Color Samples By" is a different section and stays.

**2. INFO field selection.** Selecting a field never opens a color dialog. Non-numeric fields take their scheme color, else a palette color per value (`VariantTrack.getAttributeColor`). In `defineScaleIfNeeded` (line 232, called from `getColorByInfoFieldItem`), a numeric field that no scheme, track scale or track colors cover gets a default scale applied to the track without asking: `ContinuousColorScale` from `DEFAULT_SCALE_MIN_COLOR` to `DEFAULT_SCALE_MAX_COLOR` over the loaded range (`getAttributeRange`, with `rangeEnd` for a single repeated value), set with `track.setColorScale`. Nothing is saved to preferences. Remove the `HeatmapLegendEditor` dialog, its "Update Preferences" checkbox and the discrete-scheme branch from that method, and update its javadoc. The scale is applied only if at least one loaded value can be read as a number; if none can (nothing loaded, or values that are not numbers, such as multi-valued records), do not invent a range: show a message that no loaded features have a numeric value for the field (a message box, not a color dialog) and leave the selection alone. After selection, the colors are changed from the "INFO Field <field> Colors..." item (`VariantColorLegendDialog`), whose "Edit Scale..." button edits the range and colors, as before.

**2a. "Values Are Categories..." in the legend dialog.** The first-time dialog used to offer the choice that a numeric field holds categories; with it gone, `VariantColorLegendDialog` offers it instead.
* Add a "Values Are Categories..." button beside "Edit Scale...", visible when the displayed scale is a `ContinuousColorScale` and `track.isNumericAttribute(infoKey)` (set in `populate()`). Allele frequency fields (`ColorStopScale`) do not get it.
* On click: refuse with the existing "more than 500 distinct values" message if the loaded values exceed `VariantTrackMenuHelper.MAX_CATEGORICAL_VALUES`; otherwise ask for confirmation, since the change applies at once rather than on OK and there is no way back in the dialog.
* Then give each loaded value a distinct color with `DistinctColors.next`, remove the track scale (`track.setColorScale(infoKey, null)`), set the colors with `track.setAttributeColor`, clear `pendingScale` and `pendingColors`, and repopulate, so the dialog shows a row per value to click and change.
* A scheme with a scale for the field (`VariantColorSchemes.getSchemeForScale`) would still win over the track, so save a categorical scheme (`VariantColorSchemes.saveScheme(infoKey, colors, true)`, which replaces the field's scheme) whenever "Update Preferences" is checked or such a scheme exists. If the save fails, show the error and change nothing.
* Going back from categories to a scale is not offered; removing the field's saved scheme restores the default scale.

**3. Matching check mark.** The "INFO Field" `JMenu` can't show a real check box, so it uses the hand-drawn `CheckIcon` (line 501). Make it paint the look-and-feel's own mark, `UIManager.getIcon("CheckBoxMenuItem.checkIcon")`, against a private selected `JCheckBoxMenuItem`, with the size taken from that icon and the foreground copied from the painted component. The menu is rebuilt on every popup, so fetching the icon in the constructor follows theme changes, including dark mode. If a look-and-feel has no such icon (null), or its icon paints nothing for a detached item (Nimbus's does), fall back to the current hand-drawn 12x12 mark instead. Decide this once when the `CheckIcon` is created, by painting the standard icon into a small image and checking for any painted pixel, rather than by naming look-and-feels. Check by eye on macOS in light and dark, next to the "Allele Frequency" and "Allele Fraction" marks; if the standard icon does not paint on a submenu, fall back to a custom icon reading `UIManager.getColor("CheckBoxMenuItem.foreground")`.

**4. One color item.** Replace `getAlleleFrequencyColorsItem` (line 160) and `getEditInfoColorsItem` (line 372) with one `getEditColorItem(track)`, added once after the INFO Field entry whether or not the file has INFO fields. A color-by option is always selected, so the item is never disabled. It opens:
* an INFO field selected (`ATTRIBUTE` mode with a non-null `getColorByAttribute()`): `VariantColorLegendDialog(frame, track, infoKey)`;
* otherwise, with `track.isAlleleFrequencyBars()` true: `AlleleFrequencyBarColorsDialog(frame, track)`;
* otherwise: `AlleleFrequencyColorsDialog(frame, track, colorBy)`, with `COLOR_BY_ALLELE_FRACTION` in fraction mode and `COLOR_BY_ALLELE_FREQUENCY` otherwise.

Bars vs color scale follows the "Show allele frequency / fraction as" preference through `isAlleleFrequencyBars()`. A choice saved in a session for the track still overrides it (decided earlier).

**5. Titles and item text.** Menu item text, built from the same state: `"INFO Field " + infoKey + " Colors..."` for an INFO field, otherwise `"Allele Frequency/Fraction Colors..."` (no "Edit"; the ellipsis stays). Dialog titles, set in each `super(owner, ...)` call: `AlleleFrequencyBarColorsDialog` (line 45) and `AlleleFrequencyColorsDialog` (line 42) become "Allele Frequency/Fraction Colors"; `VariantColorLegendDialog` (line 63) becomes `"INFO Field " + infoKey + " Colors"`. Menu and dialogs are rebuilt on each popup and click, so they follow the current selection with no extra code. The warning popup in `AlleleFrequencyColorsDialog` keeps its own "Allele Frequency Colors" title.

**6. Bars dialog.** `AlleleFrequencyBarColorsDialog.java` and `VariantRenderer.java`:
* Remove the "No value or 0%" example bar (line 208), make "Reference block" index 2, and reduce the preferred width from `4 * SPACING` (line 196) to `3 * SPACING`. Update the `SampleBars` javadoc.
* Renderer: in the bars branch (`percent <= 0`, line 170), a variant with no value or 0 currently gets `colorAlleleRef`, which is the track color. For a variant that is not a reference block, use `colorAlleleBandRef` (`colorAlleleBandRefAlpha` when filtered), the same reference color as other bars. A gVCF reference block (`variant.isNonRef()`) stays gray. This changes how those bars look on existing tracks.
* Replace the top text (line 83 onward; delete the `value` string) with one fixed message, `<br>` after each sentence: "The lower part of each bar is colored in the variant color and its height is in proportion to the allele frequency or fraction." / "Allele frequency = AF; or GMAF if a variant has no AF." / "Allele fraction = AC divided by AN." / "A variant with no value is drawn in the reference color." / "A gVCF reference block is drawn in gray." / "Filtered variants are drawn faded." The old "alleles are added" and "Changes apply to this track" sentences go.
* Add a test in `AlleleFrequencyColorsTest` (it already renders site bands): in bars mode the no-AF variant (`allele_frequency_sites.vcf`, index 4) paints `track.getAlleleFrequencyReferenceColor()`, and the reference block (index 5) stays `new Color(200, 200, 215)`.

**7. Color scale dialog.** `AlleleFrequencyColorsDialog.java`, line 59: replace the top text with "Variants are colored by allele frequency or fraction, blended on a log scale between stops.<br>Variants with no value take the color of the lowest stop." The old "Changes apply to this track." goes. Display only.

### Variant preferences

The tab is built by `PreferencesEditor.init` from `preferences.tab`: each `##` heading is a titled box and a `---` line is a horizontal separator inside it. The INFO colors box is not in `preferences.tab`; `PreferencesEditor` adds `VariantColorSchemePanel` after the rows, as its own box.

**1. `preferences.tab`.** Replace `##Variant Colors` (line 65) with two boxes, keeping keys, types and defaults, so saved preferences are unaffected:
* `##Genotype colors`: `HOMREF.COLOR`, `HOMVAR.COLOR`, `HETVAR.COLOR`, `NOCALL.COLOR`, unchanged.
* `##Variant site colors`: the remaining rows in their current order, `AF_REF.COLOR`, `AF_VAR.COLOR`, `VARIANT_COLOR_BY_ALLELE_FREQ`, `VARIANT_ALLELE_FREQUENCY_DISPLAY`, then a `---` line. Relabel the top two: "Allele frequency / fraction shown as bar height - Reference" and "Allele frequency / fraction shown as bar height - Variant".
* The tooltip on `VARIANT_ALLELE_FREQUENCY_DISPLAY` says "the reference and variant colors above", and they are still above it, so it stays correct. Nothing else refers to the "Variant Colors" heading.

**2. `PreferencesEditor.java`** (the `Variants` special case, lines 324-330). `group` still refers to the last box built, "Variant site colors". Add the scheme panel to that group's grid at the next `row` (width 4, like the separator, `weightx` 1, horizontal fill) and `group.add(schemePanel)`, instead of adding it as its own row of `contentGrid`. Remove the two `contentGrid`/`content` lines and update the comment. Check the table width inside the box and that the taller tab still scrolls.

**3. `VariantColorSchemePanel.java`.**
* Delete `setBorder(createTitledBorder("Variant INFO colors"))` (line 47); the enclosing box supplies the title. The panel has one caller.
* Description label (lines 49-52): remove the first line ("Colors for the values of VCF INFO attributes, used by the track menu's Color By > INFO Field option."), and in the second line change "one file per attribute" to "one file per VCF INFO attribute". The third line ("Importing a file replaces the colors for every attribute it covers.") is unchanged. The table and buttons are not touched.

### Testing

Menu, dialog and preferences changes are display only: open a VCF with INFO fields and check items 1-7, then View > Preferences > Variants (two boxes; four genotype swatches; the four site rows, a line, then the INFO table and buttons; a changed swatch and the display selection still save). The one logic change, the renderer in item 6, gets the test described there. Finish with `./gradlew compileJava` and the two test runs.

## Implementation status

All of the plan is implemented on branch `variant-info-colors`, uncommitted. It compiles, and the `org.igv.variant` and `org.igv.prefs` tests pass; the baseline run before the changes also passed. The menu, dialogs and preferences tab have not been checked in the running app.

**Done**
- **Menu items 1-5:** The label, the INFO field behavior, the check mark, the single color item and its text, and the dialog titles. The INFO field wording is "INFO Field <field> Colors...".
- **Check mark on Nimbus:** The standard check mark painted nothing under Nimbus (0 pixels, against 57-169 under Metal, Mac OS X and Motif in a headless test). `CheckIcon` now paints it once into an image when created and uses the drawn mark if nothing appears. Not checked in the app; on Nimbus the mark is the simple drawn one, so it may not match the others exactly.
- **Item 2:** Selecting a numeric INFO field applies a default scale with no dialog, only when at least one loaded value is numeric; otherwise a message is shown and the selection is unchanged.
- **Item 2a:** The "Values Are Categories..." button is in the INFO field colors dialog.
- **Items 6 and 7:** The bars dialog and the color scale dialog have the new text. The "No value or 0%" bar is gone, and `VariantRenderer.java` draws variants with no value or 0 in the Reference color.
- **Variant preferences:** "Genotype colors" and "Variant site colors" boxes, the separator, the relabeled rows, and the INFO colors table inside the site colors box with its shortened description.
- **Test added:** `testBarsNoValueIsReferenceColor` in `AlleleFrequencyColorsTest`. No other tests were added or changed.

**Known leftovers from the code review**
- **Unused code:** `saveDiscreteScheme` and the `HeatmapLegendEditor` import in `VariantTrackMenuHelper.java` are no longer used.
- **Unsaved edits:** "Values Are Categories..." discards unsaved edits in the INFO field colors dialog, and applies at once rather than on OK.
- **Categories choice:** With "Update Preferences" unchecked and no values loaded, nothing records the choice, so the default scale would come back on the next selection.

## IGNORE FOR NOW - Consider for later
* Add menu item to switch between bars and color scale for variant allele / frequency
* Add menu item to color sample genotypes by: genotype or type
* Make the menu specific to variant sites or sample genotypes, depending on where the user clicked (like in alignment tracks vs coverage track)


