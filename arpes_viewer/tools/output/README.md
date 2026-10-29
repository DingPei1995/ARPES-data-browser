# Figures, export and the experiment log: technical reference

The figure composer's model and file, export (PNG/PDF/SVG, Igor, `.nxs`), and the experiment log: `figure.py`, `figfile.py`, `export.py`, `igor.py`, `explog.py` in this folder.

## Figures and export

The **figure composer** is laid out like a manuscript page. With nothing
selected it opens a blank page of empty panels (a)–(d); each panel can be
dragged and resized on the page (snapping to the others' edges), filled
from the main list or with a picture from a file, and given an inset -- a
magnified region of itself (*Zoom-in inset*) or another panel drawn inside
it. Arrows, lines, boxes, ellipses and text are drawn straight on the page
with the tools on its left, and edited (colour, width, fill, text) when
selected; Undo/Redo covers it all. The figure tools beside the page set
the selected panel's axes, fonts, colours and the export; shared colour
scales, labels only on the outer edges and a slice-series page (a row of
constant-energy contours, or of cuts) from one cube work as before. It is
reached from the list (tool panel → Advanced → **Figure composer...**),
from a frozen slice's **Figure composer...** button, or from a map's
**Slice series figure...**. Empty slots are never written to a file.
Several panels can be selected (Ctrl/Shift-click) and lined up, spaced
or matched in size; the letters line up and run (a), (b), (c). **Save
as...** writes the composer's own `.afig` file, which reopens for editing,
and beside it a folder *data for <name>* with each panel's data as shown and
the original measurement it came from.

Export gives what the panel actually shows, in the formats that make sense
for it -- image, data, or both. Data go to a `.nxs` this program reopens, or
**for Igor Pro**: an Igor text file (`.itx`) or a binary wave (`.ibw`,
version 5), with the wave scaling, units, dimension labels and the
metadata in the wave note (`tools/output/igor.py`); an unevenly spaced axis
is written as a wave of its own beside the data. The Datasets box's
**Save** offers the same two formats for whole datasets (cuts, maps,
spatial scans, curves); an Igor export is a copy, so the rows keep pointing
where they did.

## The experiment log

Two buttons side by side in the tool panel's *Datasets* box write the
beamtime's log -- one CSV row per measurement, like the hand-kept DLS and
ANTARES logs -- from the files' own metadata.

- **Log format...** sets the log up, once:
  - **Log file** -- a folder and a file name. Pointed at a log that already
    exists, the window shows *that log's* format, to adjust.
  - **Format** -- the separator and the decimal mark. By default they are
    Excel's on this computer: Excel splits a CSV at the system's list
    separator, `;` wherever the decimal mark is a comma (much of Europe), so
    a comma file opens there with every row in the first column. On
    Windows they are read from the system settings; on macOS (Excel for
    Mac, Numbers) and Linux (LibreOffice) from the decimal mark. Numbers
    are written with the decimal mark, so a spreadsheet reads them as
    numbers rather than text or dates. Ranges are written as
    `start/stop/step` or `start:step:stop`.
  - **Columns** -- pick a dataset of the main list as the reference: its
    fields are listed with their values -- the file name, data type, start
    time and the scanned ranges worked out from the data, then every
    metadata field and motor position the file has. Tick a field to add it
    as a column; rename the columns, and put them in any order (Top / Up /
    Down / Bottom). *Add empty column* gives a column to fill in by hand
    (comments). *Save format as...* / *Load format...* keep one per
    beamline.
  - **Save format** keeps it (`~/.arpes_viewer/explog_settings.json`). If
    the log already holds records and its format changed, the log is
    rewritten under the new header, every record moved with its columns
    (followed through renames and moves); columns taken out of the format
    are named first, since what they hold goes. The header is only ever
    the first row.
- **Auto Exp. log** -- one click: the selected datasets are added to the log
  in its format, after its last record (blank rows a spreadsheet left at
  the end are dropped), skipping files already in its file-name column. The
  first click, before any format was saved, opens *Log format...*. *Keep
  recording* (in the format window) logs every file added to the list from
  then on, until the program closes: add each file as it is measured and
  the log keeps itself. A log continued on another computer keeps its own
  separator.
- **What a position column says** depends on the measurement, through each
  column's *If scanned* rule:
  - a cut: the sample stood still -- the motor's value;
  - a SPEM scan: the scanned stage axis as `start/stop/step` (*SPEM:
    scanned X / Y*, or *auto* when the axis carries the motor's name);
  - a map turned with the manipulator rather than the deflector: the angle
    as `start/stop/step` (*map: scanned angle*, or *auto* by name), and a
    sample position the beamline corrected along with the angle as
    `first..last`. Diamond files keep these moving positions (`<motor>
    range` among the motors); a value stored per point is summarised the
    same way;
  - a kz map: the photon energy range (*kz map: photon energy*).

New logs start with a byte-order mark, so that Excel reads hν and Δ. A log
open in Excel cannot be written on Windows: close it there first.

The work is in `tools/output/explog.py` (Qt-free), the window in
`ui/panels/explog.py`.
