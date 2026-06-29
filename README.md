# ephys-pipe

Streamlined electrophysiology pipeline (started 2026): align MonkeyPsych behavior with TDT ephys blocks, attach UltraSort spike times, build event-aligned rasters, and plot condition-split PSTHs.

---

## Expected inputs

Everything the pipeline needs falls into four groups. Paths below use default roots from `epp_general_settings.m` (`Y:\...` on the lab machine); override via `cfg.roots.*` in project settings.

### 1. Settings

Project and version configuration live in a **sibling** `Settings` repo (resolved relative to this repo):

```
Settings/
  <project>/
    epp/
      epp_project_settings.m      ← monkeys, datasets, WINDOWS, CONDITIONS, …
      epp_<version>_settings.m    ← optional version overrides (e.g. epp_version_a_settings.m)
```

Loaded by `epp_load_cfg(project, version)` before any stage runs. See [Configuration](#configuration).

### 2. Metadata — electrode depths

One MATLAB script per monkey, co-located with the sorting table:

```
<cfg.roots.sorting_tables>/<monkey>/Electrode_depths_<Mon>.m
```

Example: `Y:\Data\Sorting_tables\Flaffus\Electrode_depths_Fla.m`

Executed with `run()`. Must define workspace variables:

| Variable | Meaning |
|----------|---------|
| `Session` | Session date(s) — scalar or array (`YYYYMMDD`) |
| `block` | Ephys block number(s) for each depth entry |
| `channels` | Channel number(s) at that depth |
| `z` | Electrode depth(s) in µm (same length as `channels` per entry) |

Used by `epp_update_sorting_table` to assign `Blocks` and `Filenumber` per channel. `channels` and `z` must match in length per entry.

### 3. Data

Raw recordings and spike-sorting outputs on disk (assumed present before preprocessing).

**TDT ephys blocks**

```
<cfg.roots.ephys_tanks>/<monkey>_phys/<YYYYMMDD>/Block-<N>/
```

Read with `TDTbin2mat_working` (epocs: `SVal`, `Tnum`, `RunN`, `Sess`). Run numbers for behavior matching: `unique(ephys_data.epocs.RunN.data)`.

**MonkeyPsych behavior**

```
<cfg.roots.behavior>/<monkey>/<YYYYMMDD>/<Mon>YYYY-MM-DD_<RR>.mat
```

`<Mon>` = first three characters of monkey name (`Flaffus` → `Fla`). `<RR>` = zero-padded run number (`01`, `02`, …), must match ephys `RunN` for that block. Variable required: `trial` (struct array).

Example: `Y:\Data\Flaffus\20180525\Fla2018-05-25_03.mat`

**UltraSort spike files (sortcode information)**

```
<cfg.roots.ephys_tanks>/<monkey>_phys/<YYYYMMDD>/dataspikes_rb<FFF>_ch<CCC>_negthr.mat
```

`<FFF>` = file number (3 digits), `<CCC>` = channel (3 digits). Variables: `cluster_class` (col 1 = sort code, col 2 = spike time in **ms**), `par.segmentends`, `par.sr`.

### 4. Sorting table (Excel) — separate from raw data

Curated unit list; **not** generated entirely by the pipeline. Lives next to the electrode-depths script:

```
<cfg.roots.sorting_tables>/<monkey>/<Mon>_UltraSort.xlsx
```

Example: `Y:\Data\Sorting_tables\Flaffus\Fla_UltraSort.xlsx`

The pipeline reads sheet **`final_sorting`** only (`epp_sorting_table_to_units`). `epp_update_sorting_table` writes **`automatic_sorting`**; you promote that to `final_sorting` by hand after review.

**Columns filled automatically** (from spikes + electrode depths): `Session`, `Filenumber`, `Blocks`, `Channel`, `z`, `Unit`, `N_spk`, `Neuron_ID`, `Site_ID`, …

**Columns you maintain manually** — preserved across updates when session/site/unit match:

| Column | Why it matters |
|--------|------------------|
| **`Dataset`** | Numeric study/dataset ID. Filtered by `cfg.datasets` in project settings (e.g. `cfg.datasets = [3]`). Rows not in that list are skipped. **Set this before running the pipeline.** |
| **`Hemisphere`** | Recorded hemisphere per unit (`L` / `R`, case-insensitive). **Required** — `epp_build_rasters` errors if missing or invalid. Used to flip lateral trial variables (see [Trial enrichment](#trial-enrichment-stage-2)). |
| **`Perturbation`** | Perturbation condition for this unit. Copied to `trial.perturbation` on every raster trial (default `0` if empty). Use for `cfg.CONDITIONS` or downstream grouping. |
| `SNR`, stability, ranking | Quality / inclusion criteria |
| Grid-hole / anatomical fields | Any extra columns you use for curation |

**Other columns the pipeline depends on:**

| Column | Use |
|--------|-----|
| `Session` | Recording date (`YYYYMMDD`) |
| `Blocks` | Ephys blocks for this unit, e.g. `\|2\|5\|7\|` |
| `Channel`, `Filenumber`, `Unit` | Spike file lookup and sort code |
| `Neuron_ID` | Output filenames for rasters / PSTHs |

Rows with `Unit == 0` are dropped.

---

## Before running the pipeline

Assuming TDT blocks and behavioral `.mat` files are already on disk:

### 1. Create the electrode depths file

Write `Electrode_depths_<Mon>.m` for each monkey (see [Metadata](#2-metadata--electrode-depths)). One script per monkey; can cover many sessions/blocks.

### 2. Run UltraSort

Sort spikes and export `dataspikes_rb<FFF>_ch<CCC>_negthr.mat` into the session folder under `<monkey>_phys/<YYYYMMDD>/`.

### 3. Update the UltraSort table

```matlab
epp_update_sorting_table('Flaffus');                          % all sessions in source folder
epp_update_sorting_table('Flaffus', [20160608, 20160609]);    % specific dates only
```

This step is part of the ephys-pipe workflow but **not** called by `epp_initiation`.

**What it does:**

1. Reads existing `<Mon>_UltraSort.xlsx` (prefers `final_sorting`, else `automatic_sorting`)
2. Loads `Electrode_depths_<Mon>.m`
3. Scans UltraSort spike `.mat` files per session; assigns `Filenumber` from depth steps (>50 µm → new file), `Blocks` from electrode depths, sort codes from `cluster_class`
4. Merges with existing manual columns (SNR, stability, **`Dataset`**, etc.) where session/site/unit match
5. Writes sheet **`automatic_sorting`** (never overwrites `final_sorting` directly)

**Spike source for updates:** merge-data folder on the lab server (`.../spikesorting/testdata/merge_data_structure/TDTbrain/<YYYYMMDD>/dataspikes*negthr.mat`). Align that path with your tanks or adjust `epp_update_sorting_table` if your layout differs.

### 4. Review and promote → ready to go

1. Open `<Mon>_UltraSort.xlsx`, compare `automatic_sorting` to `final_sorting`
2. Fill or correct **manual columns**, especially **`Dataset`**, **`Hemisphere`**, and **`Perturbation`**
3. Copy/replace the `final_sorting` sheet from `automatic_sorting` (or merge selectively)
4. Confirm `cfg.datasets` and `cfg.monkeys` in project settings match the rows you kept

You are ready when `final_sorting` reflects the units you want and paths in [Expected inputs](#expected-inputs) resolve on your machine.

---

## Running the pipeline

Add this repository to your MATLAB path, then:

```matlab
epp_initiation('project_tdt_brain', 'version_a');
```

`project` and `version` select settings under `Settings/<project>/epp/`. Stages can also be run individually:

```matlab
cfg = epp_load_cfg('project_tdt_brain', 'version_a');
epp_prepare_blocks(cfg);
epp_build_rasters(cfg);
epp_compute_statistics(cfg);
epp_plot_psth(cfg);
epp_plot_population_psth(cfg);
```

### What `epp_initiation` does

Loads `cfg` once via `epp_load_cfg`, then runs five stages:

#### Stage 1 — `epp_prepare_blocks`

1. Read `final_sorting` → build `unit_info`
2. For each relevant `(session, block)` from unit rows:
   - Load TDT ephys block
   - For each run in `epocs.RunN`: load matching behavior `.mat`, enrich trials (`MP_add_saccades_and_reaches`), synchronize (`epp_synchronization`), concatenate runs
   - Save one `block_payload` per block
3. Save `unit_info.mat` and `prepare_blocks_report.txt`

Behavior is mapped into ephys block time using state 2 as the trial anchor.

#### Stage 2 — `epp_build_rasters`

1. Load `unit_info.mat` and saved block payloads
2. Per unit: load spike file, offset behavioral timestamps by block segment start (`par.segmentends`), concatenate across blocks
3. **`epp_enrich_trials_for_unit`** — attach `perturbation`, flip lateral fields by `Hemisphere` (see below)
4. Per `cfg.WINDOWS` entry: align spikes to state events, bin, attach trial metadata
5. Save `<Neuron_ID>_raster.mat`

#### Trial enrichment (stage 2)

`epp_enrich_trials_for_unit(trials, unit_row)` runs once per unit after blocks are concatenated. Sorting-table columns (exact header names):

| Excel column | Trial field | Rule |
|--------------|-------------|------|
| `Perturbation` | `perturbation` | Copied to every trial; default `0` if empty |
| `Hemisphere` | — | `L`/`l` → multiplier **+1**, `R`/`r` → **−1**; **error** if neither |

**Lateral flip** (multiply by hemisphere multiplier):

- **Target positions and `hemifield`** (`tar_pos`, `nct_pos`, `fix_pos`, `reach_tar_pos_closest`, `saccade_tar_pos_closest`, `hemifield`): `real` part `*= multiplier`, `imag` unchanged; real-only fields stay real
- **Hands** (`demanded_hand`, `used_hand`, `reach_hand` if present): `1 → −1`, `2 → +1`, then `*= multiplier`

#### Stage 3 — `epp_compute_statistics`

1. Per unit: load aligned trials + spikes (same path as stage 2), enrich trials
2. Per trial × epoch: mean firing rate (Hz) in epoch window relative to align state
3. Per condition:
   - **Baseline comparisons** — paired test (`cfg.statistics.baseline_test`, default `paired_ttest`) of each epoch vs its `baseline` epoch
   - **Condition comparisons** — unpaired test (`cfg.statistics.condition_test`, default `unpaired_ttest`) for all condition pairs, per epoch
4. Save `unit_statistics.mat` (N-unit struct array, one field per comparison) and `unit_statistics.xlsx` under `cfg.roots.statistics`

Comparison field names: `<Condition>_<Epoch>_vs_baseline_<Baseline>` and `<Epoch>_<CondA>_vs_<CondB>` (valid MATLAB identifiers).

#### Stage 4 — `epp_plot_psth`

1. Load all raster files
2. For each window, split aligned events (trials) by `cfg.CONDITIONS` — see [Condition parameters](#condition-parameters-cfgconditions)
3. Smooth raster rates, rebin to PSTH bins, plot PSTH + condition-colored raster
4. Save `<Neuron_ID>_psth.png` and `<Neuron_ID>_psth.mat`

#### Stage 5 — `epp_plot_population_psth`

1. Load all unit `*_psth.mat` files from stage 4
2. Per window and condition: mean across units of per-unit PSTH curves; SEM across units
3. Plot PSTH only (no raster), shaded ± SEM
4. Save `population_psth.png` and `population_psth.mat` under `cfg.roots.population_psth`

### Outputs

Written under `cfg.roots.project_version` (default `Y:\Projects\<project>\<version>/`):

| Stage | Folder (`cfg.roots.*`) | Main files |
|-------|------------------------|------------|
| 1 | `processed_trials` (`.../behavior/`) | `<monkey>_<YYYYMMDD>_Block-<NNN>.mat`, `unit_info.mat`, `prepare_blocks_report.txt` |
| 2 | `raster` (`.../unit_rasters/`) | `<Neuron_ID>_raster.mat` |
| 3 | `statistics` (`.../statistics/`) | `unit_statistics.mat`, `unit_statistics.xlsx` |
| 4 | `psth` (`.../unit_psth/`) | `<Neuron_ID>_psth.png`, `<Neuron_ID>_psth.mat` |
| 5 | `population_psth` (`.../population_psth/`) | `population_psth.png`, `population_psth.mat` |

`prepare_blocks_report.txt` lines: `ephys_block - behavior_file - run=N - messages` (sync reports, missing files, empty trials, etc.).

---

## Configuration

`epp_load_cfg(project, version)` builds `cfg` in this order:

| Step | File | Required |
|------|------|----------|
| 1 | `ephys-pipe/epp_general_settings.m` | yes |
| 2 | `../Settings/<project>/epp/epp_project_settings.m` | warned if missing |
| 3 | `../Settings/<project>/epp/epp_<version>_settings.m` | warned if missing |

Each script is `run()` in the MATLAB workspace and assigns fields to `cfg` (override only what you need). `epp_load_cfg` creates any missing output root folders.

### `cfg.roots`

| Field | Default role |
|-------|----------------|
| `settings` | `../Settings/<project>/epp/` |
| `sorting_tables` | UltraSort Excel + electrode-depth scripts |
| `ephys_tanks` | TDT tanks + spike `.mat` files |
| `behavior` | MonkeyPsych `.mat` files |
| `project_version` | `Y:\Projects\<project>\<version>/` |
| `processed_trials` | `.../behavior/` |
| `raster` | `.../unit_rasters/` |
| `psth` | `.../unit_psth/` |
| `population_psth` | `.../population_psth/` |
| `statistics` | `.../statistics/` |

### Analysis parameters (typical project overrides)

- `cfg.monkeys` — e.g. `{'TDTbrain'}`
- `cfg.datasets` — numeric filter on sorting-table **`Dataset`** column
- `cfg.WINDOWS` — struct array: `name`, `align_state`, `t_start_s`, `t_end_s`
- `cfg.CONDITIONS` — PSTH / raster grouping by behavioral trial fields (see below)
- `cfg.raster_bin_size_s` — raster bin width (default 1 ms)
- `cfg.psth.bin_size_s`, `cfg.psth.smoothing_kernel`, `cfg.psth.smoothing_width_s`
- `cfg.EPOCHS` — cell table in settings (`name | state | start | end | baseline`); struct array after `epp_load_cfg`
- `cfg.statistics.baseline_test` — default `paired_ttest`
- `cfg.statistics.condition_test` — default `unpaired_ttest`

### Condition parameters (`cfg.CONDITIONS`)

Used in stage 4 (`epp_plot_psth`). Each aligned event carries the behavioral `trial` struct active at align time (from stage 2). Conditions select which events go into each PSTH curve and raster stripe color.

**Struct fields per condition:**

| Field | Meaning |
|-------|---------|
| `name` | Label in legend and saved `psth_data` |
| `color` | RGB line / raster-dot color (`0–1` or `0–255`) |
| `parameters` | Struct of trial-field filters (see matching rules below) |

**Example** (left vs right target, from `epp_project_settings.m`):

```matlab
cfg.CONDITIONS = struct( ...
    'name', {'Left', 'Right'}, ...
    'color', {[200 55 12], [100 55 22]}, ...
    'parameters', num2cell(struct('choice', {0, 0}, 'hemifield', {-1, 1})));
```

Here `Left` keeps events whose trial has `choice == 0` **and** `hemifield == -1`; `Right` uses `choice == 0` **and** `hemifield == 1`.

#### Matching rules

For each condition, `epp_plot_psth` builds a boolean mask over aligned events:

1. **Field names** in `parameters` must match fields on the behavioral `trial` struct (e.g. `choice`, `hemifield`, `target`, `correct`, … — whatever MonkeyPsych wrote into `trial`).
2. **Within one field** — the value in settings is a scalar or vector. An event matches if its trial value is **any** of the listed values (`ismember`).  
   Example: `'choice', [0 1]` pools correct and error trials into one condition.
3. **Across fields** — all listed fields must match (**AND**).  
   Example: `'choice', 0, 'hemifield', [-1 1]` → `choice == 0` and hemifield is either −1 or +1.
4. **No fixed schema** — add or omit fields freely; only the fields you list are checked. Omitted trial variables are not filtered on.

So conditions are flexible: you define arbitrary trial-field combinations, and each field accepts one or many allowed values without writing separate code paths.

#### Overlap between conditions

- **PSTH curves** — masks are computed independently. The same event can contribute to **multiple** conditions if it satisfies more than one definition (overlap is allowed).
- **Raster panel** — each event is colored by the **first** condition in `cfg.CONDITIONS` order that matches; events matching none are drawn in gray.

`n_trials` saved per condition is the number of events passing that condition’s mask (after overlap, counts need not sum to total events).

---

## Main functions

| Function | Role |
|----------|------|
| `epp_update_sorting_table` | Preprocessing: regenerate `automatic_sorting` from spikes + electrode depths |
| `epp_initiation` | Run all five pipeline stages |
| `epp_load_cfg` | Build `cfg` from settings hierarchy |
| `epp_prepare_blocks` | Stage 1: sync behavior ↔ ephys, save blocks |
| `epp_build_rasters` | Stage 2: per-unit event-aligned rasters |
| `epp_enrich_trials_for_unit` | Stage 2: copy `Perturbation`, flip lateral trials by `Hemisphere` |
| `epp_compute_statistics` | Stage 3: epoch firing-rate statistics per unit |
| `epp_plot_psth` | Stage 4: condition-split PSTHs |
| `epp_plot_population_psth` | Stage 5: population PSTH across units |
| `epp_synchronization` | Align one behavioral run to one ephys block |
| `MP_add_saccades_and_reaches` | Enrich trials with movement timing and targets |
| `epp_sorting_table_to_units` | Parse `final_sorting` into unit structs |
