"""The repo's pointers to the config builder stay true."""

import re

from gui_support import GUI_PATH, REPO_ROOT

PAGES_ROOT = "https://hawkjo.github.io/bcwithqc/"


def test_readme_links_to_the_config_builder_on_github_pages():
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    links = re.findall(r"\[Open the config builder\]\((.*?)\)", readme)
    assert links, "README.md has no 'Open the config builder' link"
    for link in links:
        assert link.startswith(PAGES_ROOT)
        assert (REPO_ROOT / link[len(PAGES_ROOT) :]) == GUI_PATH
    assert GUI_PATH.exists()


def test_github_pages_serves_files_as_they_are():
    """An empty .nojekyll at the repo root stops GitHub Pages from running Jekyll."""
    assert (REPO_ROOT / ".nojekyll").exists()
