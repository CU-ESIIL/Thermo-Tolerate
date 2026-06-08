# Thermo ToleRate

**Redefining ectotherm thermal tolerance in a variable climate using big data**

An ESIIL Working Group building a time-dependent framework for ectotherm thermal tolerance
that accounts for both temperature intensity and exposure duration.

This repository is one connected system:

- The **repository** is where the science happens — data, notebooks, workflows, and outputs.
- The **website** (built from `docs/`) is where the science is shared.
- **GitHub** connects them through commits and version history.

## Goals

1. Build a harmonized, open database of ectotherm thermal tolerance experiments with
   duration- and performance-based outcomes.
2. Generate thermal tolerance surfaces from physiological data and link them to climate records.
3. Test whether these surfaces predict heatwave-driven population mortality across
   terrestrial and aquatic taxa.

## Repository layout

```text
.
├── README.md       # This file
├── mkdocs.yml      # Website navigation, theme, plugins, and edit links
├── docs/           # Markdown source for the public website (docs/index.md is the homepage)
├── scripts/        # Build helpers and site health checks
├── templates/      # Reusable meeting-note templates
├── containers/     # Optional runtime and environment setup
└── ...             # Add data, notebooks, workflows, and figures as the science grows
```

## Preview the website locally

```bash
pip install -r requirements.txt
python scripts/generate_image_slots.py
mkdocs serve
```

## Status

Pre-meeting coordination — preparing for Meeting 1 (Boulder, June 2026).

Project Leader: Kelsey Lyberger (kelsey.lyberger@asu.edu)
