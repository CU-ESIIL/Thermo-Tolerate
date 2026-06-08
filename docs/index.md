# Thermo ToleRate

**Redefining ectotherm thermal tolerance in a variable climate using big data**

An ESIIL Working Group building a time-dependent framework for ectotherm thermal
tolerance — one that accounts for both *how hot* it gets and *how long* the heat lasts.

- **Repository** (where the science happens — data, code, workflows, outputs): [add link]
- **Website** (where the science is shared): this site

## Abstract

Thermal tolerance determines how ectothermic animals survive climate extremes, but
traditional metrics like critical thermal maximum (CTmax) and lethal temperature (LT50)
ignore that tolerance depends on both the *intensity* and *duration* of heat exposure.
This working group redefines thermal tolerance as a time-dependent surface that integrates
temperature and exposure time, linking experimental physiology with real-world climate
variability.

We will (1) assemble a harmonized, open database of thermal tolerance experiments with
duration- and performance-based outcomes, (2) develop Bayesian hierarchical models to
estimate continuous tolerance surfaces and compare them across taxa, and (3) simulate
expected survival under realistic climate records (ERA5-Land, OISST) and test predictions
against observed abundance in NEON terrestrial insects and marine copepods.

## Research Question

How can we quantify and model thermal tolerance as a time-dependent surface that predicts
ectotherm population responses to real-world heatwave dynamics across taxa and environments?

## Goals

1. Build a harmonized, open database of ectotherm thermal tolerance experiments that
   include duration- and performance-based outcomes.
2. Develop pipelines that generate thermal tolerance surfaces from physiological data and
   link them to temporally explicit climate data.
3. Test whether these surfaces predict population mortality during heatwaves across
   terrestrial and aquatic ectotherm taxa.

## Planned Outputs

1. **Database + data paper** — a cross-taxon database of heat-tolerance records combining
   temperature intensity and exposure duration.
2. **Conceptual paper** — general patterns in the shape and scaling of time-dependent
   survival, and a shared language for comparing tolerance landscapes.
3. **Reproducible Bayesian workflow** — an R pipeline for fitting tolerance surfaces,
   archived with a DOI.
4. **Proof-of-concept paper** — applying the workflow to NEON mosquitoes and beetles and
   to marine copepod records.

## Core Taxa & Data

- **Mosquitoes & beetles** — NEON pitfall and CO₂-trap abundance, forced with hourly
  ERA5-Land temperature.
- **Copepods** — Helgoland Roads and Western Channel Observatory time series, forced with
  daily OISST.
- **Tolerance compilations** — GlobTherm, ThermoFresh, Rezende (2014), plus standardized
  records from individual studies.

## Current Phase

Pre-meeting coordination (Mar–May 2026): establishing communication, compiling datasets,
and building metadata templates ahead of **Meeting 1 (Boulder, June 2026)**.

## Team

| Name | Role | Institution |
| --- | --- | --- |
| Kelsey Lyberger | Project Leader | Arizona State University |
| _Add collaborator_ | | |

## Key Links

- Main working document: [link]
- GitHub repository: [link]
- Data / resources: [link]
- Outputs / dashboard: [link]
