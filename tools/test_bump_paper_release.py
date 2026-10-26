"""Unit checks for the release-tag URL rewrite in bump_paper_release.py."""

from __future__ import annotations

import importlib.util
import pathlib


def _load():
    path = pathlib.Path(__file__).resolve().parent / "bump_paper_release.py"
    spec = importlib.util.spec_from_file_location("bump_paper_release", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rewrite_release_urls_changes_blob_and_tree_tags_only() -> None:
    module = _load()
    text = (
        "\\url{https://github.com/FloatingPragma/observer-patch-holography/blob/r2052/docs/BASIS.md}\n"
        "\\url{https://github.com/FloatingPragma/observer-patch-holography/tree/r2052/model}\n"
        "release r2052 was deposited; https://github.com/FloatingPragma/observer-patch-holography/releases/tag/r2052\n"
        "\\url{https://github.com/FloatingPragma/observer-patch-holography/blob/r2051/README.md}\n"
    )
    out = module.rewrite_release_urls(text, "r2052", "r2053")
    assert "blob/r2053/docs/BASIS.md" in out
    assert "tree/r2053/model" in out
    assert "release r2052 was deposited" in out
    assert "releases/tag/r2052" in out
    assert "blob/r2051/README.md" in out
    assert "r2052/" not in out.replace("releases/tag/r2052", "")


def test_plan_release_url_rewrites_skips_archive(tmp_path: pathlib.Path) -> None:
    module = _load()
    (tmp_path / "paper" / "archive").mkdir(parents=True)
    live = tmp_path / "paper" / "live.tex"
    live.write_text("blob/r2052/a blob/r2052/b", encoding="utf-8")
    live.write_text(
        "https://github.com/FloatingPragma/observer-patch-holography/blob/r2052/a "
        "https://github.com/FloatingPragma/observer-patch-holography/blob/r2052/b",
        encoding="utf-8",
    )
    (tmp_path / "paper" / "archive" / "old.tex").write_text(
        "https://github.com/FloatingPragma/observer-patch-holography/blob/r2052/c", encoding="utf-8"
    )
    plan = module.plan_release_url_rewrites(tmp_path, "r2052", "r2053")
    assert plan == {live: 2}
