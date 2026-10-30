"""The README on PyPI: its relative links and images made absolute."""
from devtools.pypi_readme import absolute


def test_relative_links_and_images_point_at_the_tag():
    text = ('See [the guide](docs/USER_GUIDE_en.md#2-x), ![a](docs/images/a.png),\n'
            '<img src="docs/images/b.png">, [web](https://x.org), [here](#bug-reports),\n'
            '[mail](mailto:a@b.c).')
    out = absolute(text, "v1.2.1")
    repo = "https://github.com/DingPei1995/ARPES-viewer"
    raw = "https://raw.githubusercontent.com/DingPei1995/ARPES-viewer"
    assert f"({repo}/blob/v1.2.1/docs/USER_GUIDE_en.md#2-x)" in out
    assert f"({raw}/v1.2.1/docs/images/a.png)" in out
    assert f'src="{raw}/v1.2.1/docs/images/b.png"' in out
    assert "(https://x.org)" in out and "(#bug-reports)" in out
    assert "(mailto:a@b.c)" in out
