# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

IGV (Integrative Genomics Viewer) is a desktop Java application for visualizing genomic data. It supports Mac, Windows, and Linux. Requires Java 21.

## Build & Run Commands

```bash
# Build distribution (output: build/IGV-dist/)
./gradlew createDist

# Run all tests (headless, 2GB max heap)
./gradlew test

# Run a single test class or method
./gradlew test --tests org.igv.ucsc.TrixTest
./gradlew test --tests org.igv.ucsc.TrixTest.testMethod

# Build jar only
./gradlew jar

# Platform-specific distributions
./gradlew createMacDistZip
./gradlew createLinuxDistZip
./gradlew createWinDist
```

After `createDist`, launch with `build/IGV-dist/igv.sh` (Linux), `igv.command` (Mac), or `igv.bat` (Windows).

## Test Notes

- Tests run headless (`java.awt.headless=true`); long-running tests excluded by default
- Large test data (~1GB) available separately from ftp://ftp.broadinstitute.org/pub/igv/largedata.zip — unzip to `test/largedata/`

## Architecture

**Entry point:** `org.igv.ui.Main` → creates the Swing JFrame and calls `open()` which initializes the singleton `org.igv.ui.IGV`.

**Core singleton:** `org.igv.ui.IGV` manages the main window, genome, data model, and track list. Only one IGV instance is allowed per JVM.

**Key subsystems:**

| Package | Responsibility |
|---|---|
| `org.igv.ui` | Swing GUI: main window, panels, dialogs, menus |
| `org.igv.track` | Track abstraction — all data types extend `Track`; `TrackMenuUtils` builds context menus |
| `org.igv.renderer` | Rendering engine; each track type has a corresponding renderer |
| `org.igv.feature.genome` | Reference genome management (`Genome`, `GenomeManager`) |
| `org.igv.session` | Session save/load (XML and JSON formats) |
| `org.igv.data` | Data loading and caching layer |
| `org.igv.sam` | SAM/BAM alignment handling |
| `org.igv.variant` | VCF/variant track (`VariantTrack`) |
| `org.igv.seg` | Segmentation data (`SegTrack`) |
| `org.igv.ucsc` | UCSC formats: BigBed, BigWig, TwoBit, TDF |
| `org.igv.batch` | Headless batch mode command execution |
| `org.igv.event` | Pub/sub event bus (`IGVEventBus`, `IGVEventObserver`) |
| `org.igv.prefs` | User preferences (`PreferencesManager`) |
| `org.igv.util` | HTTP, file I/O, stream utilities |
| `org.igv.aws` | AWS S3 / Cognito integration |

**UI panel hierarchy:** `MainPanel` contains a `DataPanel` (the scrollable genome view) composed of `TrackPanel` rows. Each `TrackPanel` has a `TrackNamePanel` (left label area) and a `DataPanelContainer` (right drawing area). Panels communicate via `IGVEventBus`.

**Track rendering:** `Track.render(Graphics2D g, RenderContext context)` is the core drawing method. `RenderContext` carries the current viewport locus, scale, and panel dimensions.

**Session files:** Supported in both XML (legacy `.xml`) and JSON (`.json`) formats. The `session/` package handles serialization; example sessions are in `test/sessions/`.
