"""Make README.md's relative links absolute, for the page on PyPI.

PyPI shows README.md as the project's description, but it has no repository
to resolve ``docs/...`` against: images would be missing and links dead.
The publish workflow runs this on its own copy just before building, so the
links point at the released tag's files; the README in the repository keeps
its relative links. Usage: ``python devtools/pypi_readme.py README.md v1.2.3``.
"""
import re
import sys

REPOSITORY = "https://github.com/DingPei1995/ARPES-viewer"


def absolute(text: str, ref: str) -> str:
    """``text`` with every relative link and image made absolute at ``ref``."""
    def target(path: str, image: bool) -> str:
        if re.match(r"[a-z]+:|#", path):         # http:, mailto:, an anchor
            return path
        if image:
            return f"https://raw.githubusercontent.com/DingPei1995/ARPES-viewer/{ref}/{path}"
        return f"{REPOSITORY}/blob/{ref}/{path}"

    text = re.sub(r"(!?)\[([^\]]*)\]\(([^)\s]+)\)",
                  lambda m: f"{m[1]}[{m[2]}]({target(m[3], bool(m[1]))})", text)
    return re.sub(r'src="([^"]+)"', lambda m: f'src="{target(m[1], True)}"', text)


if __name__ == "__main__":
    path, ref = sys.argv[1], sys.argv[2]
    with open(path, encoding="utf-8") as f:
        text = f.read()
    with open(path, "w", encoding="utf-8") as f:
        f.write(absolute(text, ref))
