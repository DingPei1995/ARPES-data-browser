# agentic/ — local agentic workflow for the ARPES viewer

Lets an AI agent (e.g. Claude, through the desktop app linked to your computer) process data **inside the viewer on your own machine**,
without uploading raw data or code: the agent drops small Python command files into a folder, the live viewer executes them in its
GUI process, and results come back as text logs, screenshots and figures that you can check.

Read the full guide first: [`docs/AGENTIC_WORKFLOW_en.md`](../docs/AGENTIC_WORKFLOW_en.md).

Quick start (Windows):
1. `copy config.example.json config.json` and edit `data_dir`.
2. `python launch_viewer.py` — the viewer opens (monitor 2 by default) with the data listed.
3. Optional headless jobs: double-click `agent_runner.pyw`.
4. Put command files into `viewer_cmd/` (live viewer) or `queue/` (headless); read `out/*.log`; screenshots are in `shots/`.

`recipes/` holds building blocks adapted from a real session (normalise, cut k conversion + Fermi fit + offset, arbitrary cut from a map,
two-peak MDC fit). They call the viewer's internal functions of v1.2.3 and may need changes after updates; check each on one dataset first.
Working folders are git-ignored.
