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

# Project Constraints & Rules

## Codebase Status
- This project is mature, highly stable, and has been in production for many years.
- DO NOT rewrite, refactor, or optimize any working code unless required or explicitly requested.

## Coding Style & Principles
- **Minimalist Changes:** Always choose the path that modifies the fewest lines of code necessary to accomplish the objective.
- **Scope Restriction:** Confine your modifications strictly to the files directly responsible for the task.
- **Maintain Layout & Formatting:** Respect all existing architecture, naming conventions, and code formatting rules exactly. Avoid modernizing patterns if they deviate from the legacy standard.  Making suggestions for modernization is encouraged, but do not implement them unless explicitly requested.

## Verification
- Run existing tests before making any changes to establish a baseline.
- Run tests after changes to ensure absolute zero regression.

### Testing Rules & Constraints
- **Do Not Default to Writing Tests:** Not every change or implementation requires a new unit test. Only add or update tests if introducing new business logic, complex state mutations, or fixing a bug with reproducible edge cases.
- **Banned Test Patterns:**
    - **No Trivial Tests:** Never write tests for simple getters/setters, configuration files, logger statements, or basic wrapper functions.
    - **No Tautological Tests:** Avoid tests that simply mirror the implementation code or check if a function returns what it was hardcoded to return.
    - **No "Mock-Everything" Tests:** Do not write tests where every single dependency is mocked out to the point where the test only verifies that dependencies were called.
- **Verification Rule:** If a code change does not warrant a new test, explicitly state: *"No new unit tests required because [reason]."* Do not generate placeholder, empty, or redundant test files.
  

## Test Notes

- Tests run headless (`java.awt.headless=true`); long-running tests excluded by default

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
| `org.igv.alignment` | SAM/BAM alignment handling |
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

**Track hubs:** The "Track Hubs" menu is rebuilt by `IGVMenuBar.updateHubsMenu` on every `GenomeChangeEvent`, from three sources: hubs declared in the genome json `hubs` array (loaded by the `Genome` constructor; the first is the "genome hub" supplying default annotations), the built-in ENCODE / 4DN choosers, and hubs the user selected from the UCSC registry or added by URL (`HubRegistry`, persisted in `~/igv/hubs.txt`). `HubParser.loadHubs` fetches genome hubs in parallel under a time budget shared by all of a genome's hubs (so the wait does not scale with their number) so an unresponsive hub server cannot block genome loading; hubs that fail or time out are logged with their URL and appear as disabled "(failed to load)" menu items rather than silently vanishing.

**Session files:** Supported in both XML (legacy `.xml`) and JSON (`.json`) formats. The `session/` package handles serialization; example sessions are in `test/sessions/`.

## Git commit messages

- Keep commit messages extremely concise.
- A single short subject line (under 50 characters) whenever possible; no body.
- No bulleted lists, no recaps of what changed file-by-file, no explanation of the reasoning.
- Add a body only when the *why* is genuinely non-obvious from the diff — then one or two sentences, not a summary of the change.

