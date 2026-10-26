"""Build the PDF editions of the user guides in ``docs/pdf/``.

Each ``docs/USER_GUIDE_<lang>.md`` is turned into HTML by pandoc (GitHub
flavour, so the section links resolve exactly as they do on GitHub), styled
for A4 print, and printed to PDF by headless Chromium through Playwright,
with page numbers and a bookmark outline. Links between the guides, and to
the README and the changelog, point at GitHub.

Nothing here is part of the program; run it by hand after editing a guide
or re-taking the screenshots (``devtools/guide_screenshots.py``).

Needs:  pandoc, and ``pip install playwright`` with a Chromium it can use
        (``playwright install chromium``, or set ``CHROMIUM`` to a Chrome or
        Chromium executable).

Usage:  python devtools/build_guide_pdfs.py [lang ...]     (default: all)
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(HERE, "docs")
OUT = os.path.join(DOCS, "pdf")
REPO = "https://github.com/DingPei1995/ARPES-data-browser/blob/main"

GUIDES = {
    "en": {"title": "ARPES viewer — User Guide", "lang": "en"},
    "zh": {"title": "ARPES viewer 用户指南", "lang": "zh-CN"},
    "fr": {"title": "ARPES viewer — Guide de l'utilisateur", "lang": "fr"},
}

CSS = """
@page { size: A4; margin: 16mm 15mm 18mm 15mm; }
html { font-family: "DejaVu Sans", "Noto Sans", "Helvetica Neue", Arial,
       "WenQuanYi Zen Hei", "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
       font-size: 9.6pt; line-height: 1.5; color: #1d2330; }
body { margin: 0; }
h1 { font-size: 23pt; line-height: 1.2; margin: 0 0 4mm; color: #16365f;
     border-bottom: 2.5px solid #2f6db3; padding-bottom: 3mm; }
h2 { font-size: 14.5pt; color: #16365f; margin: 7mm 0 2.5mm;
     border-bottom: 1px solid #c9d6e8; padding-bottom: 1mm;
     break-after: avoid; }
h3 { font-size: 11.5pt; color: #2f6db3; margin: 5mm 0 2mm; break-after: avoid; }
p, li { orphans: 3; widows: 3; }
a { color: #2f6db3; text-decoration: none; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 8.6pt;
       background: #eef2f8; padding: 0.5px 3px; border-radius: 2px; }
pre { background: #f3f6fb; border: 1px solid #d9e2ef; border-radius: 3px;
      padding: 2mm 3mm; break-inside: avoid; white-space: pre-wrap; }
pre code { background: none; padding: 0; }
hr { border: 0; border-top: 1px solid #d9e2ef; margin: 5mm 0; }
img { display: block; max-width: 100%; max-height: 118mm; margin: 2mm auto;
      border: 1px solid #d0d7e2; break-inside: avoid; }
p:has(> img:only-child) { break-inside: avoid; text-align: center; }
table { border-collapse: collapse; width: 100%; margin: 2.5mm 0;
        font-size: 8.8pt; break-inside: auto; }
th, td { border: 1px solid #cdd6e3; padding: 1.2mm 2mm; vertical-align: top;
         text-align: left; }
th { background: #e8eef7; }
tr { break-inside: avoid; }
td img { max-height: 95mm; margin: 0 auto; }
/* tall, narrow dialogs sit beside the text instead of taking a page each */
img.portrait { float: right; max-width: 46%; max-height: 135mm;
               margin: 1mm 0 3mm 5mm; }
td img.portrait { float: none; max-width: 100%; margin: 0 auto; }
h2, h3, hr, table, pre { clear: both; }
/* the table of contents straight after the title */
h2#contents + ol, h2#目录 + ol, h2#sommaire + ol { columns: 2; column-gap: 8mm; }
"""

FOOTER = ("<div style='width:100%;font-size:7.5px;color:#7a8597;"
          "font-family:DejaVu Sans,sans-serif;padding:0 15mm;display:flex;"
          "justify-content:space-between'><span>{title}</span>"
          "<span><span class='pageNumber'></span> / "
          "<span class='totalPages'></span></span></div>")


def _png_size(path: str):
    """(width, height) of a PNG, from its header."""
    import struct
    with open(path, "rb") as handle:
        head = handle.read(24)
    return struct.unpack(">II", head[16:24])


def _mark_portraits(body: str) -> str:
    """Give every screenshot taller than it is wide the ``portrait`` class."""
    def mark(match):
        tag, src = match.group(0), match.group(1)
        try:
            width, height = _png_size(os.path.join(DOCS, src))
        except (OSError, ValueError):
            return tag
        return tag.replace("<img ", '<img class="portrait" ', 1) \
            if height > 1.15 * width else tag
    return re.sub(r'<img [^>]*src="(images/[^"]+)"[^>]*>', mark, body)


def to_html(lang: str) -> str:
    source = os.path.join(DOCS, f"USER_GUIDE_{lang}.md")
    body = subprocess.run(["pandoc", "-f", "gfm", "-t", "html5", source],
                          check=True, capture_output=True, text=True).stdout
    # The guides link to each other and to files beside them in the repo.
    # A relative link would be printed as this machine's absolute path, so
    # the other guides are linked where they live on GitHub.
    body = re.sub(r'href="USER_GUIDE_(\w+)\.md"',
                  rf'href="{REPO}/docs/pdf/ARPES_viewer_user_guide_\1.pdf"', body)
    body = re.sub(r'href="\.\./([\w.]+\.md)"', rf'href="{REPO}/\1"', body)
    body = _mark_portraits(body)
    body = body.replace('src="images/', 'src="../images/')
    meta = GUIDES[lang]
    return (f"<!doctype html><html lang='{meta['lang']}'><head><meta charset='utf-8'>"
            f"<title>{meta['title']}</title><style>{CSS}</style></head>"
            f"<body>{body}</body></html>")


def build(langs):
    from playwright.sync_api import sync_playwright

    os.makedirs(OUT, exist_ok=True)
    launch = {}
    if os.environ.get("CHROMIUM"):
        launch["executable_path"] = os.environ["CHROMIUM"]
    written = []
    with sync_playwright() as p:
        browser = p.chromium.launch(**launch)
        page = browser.new_page()
        for lang in langs:
            # Written beside the PDFs so the relative image paths resolve.
            html_path = os.path.join(OUT, f"_build_{lang}.html")
            with open(html_path, "w", encoding="utf-8") as handle:
                handle.write(to_html(lang))
            try:
                page.goto("file://" + html_path, wait_until="networkidle")
                pdf = os.path.join(OUT, f"ARPES_viewer_user_guide_{lang}.pdf")
                page.pdf(path=pdf, format="A4", print_background=True,
                         prefer_css_page_size=True, outline=True, tagged=True,
                         display_header_footer=True, header_template="<div></div>",
                         footer_template=FOOTER.format(title=GUIDES[lang]["title"]))
                written.append(pdf)
            finally:
                os.remove(html_path)
        browser.close()
    return written


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(GUIDES)
    for path in build(wanted):
        print(path)
