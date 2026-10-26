#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path
import json


RELEASE_INFO_RELATIVE = Path("paper/release_info.tex")
CLAIM_REGISTRY_RELATIVE = Path("claims/claim_registry.yaml")
RELEASE_ID_MACRO = "OPHPaperReleaseID"
RELEASE_DATE_MACRO = "OPHPaperReleaseDate"
RELEASE_NUMBER_PATTERN = re.compile(r"^(?P<prefix>.*?)(?P<number>\d+)$")
TEX_SOURCE_DIRS = ("paper", "flagship", "extra", "cosmology", "book")
REPOSITORY_URL_FRAGMENT = "FloatingPragma/observer-patch-holography/"
NEXT_STEP = (
    "Next: run python3 tools/refresh_paper_release.py --publication"
)


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    release_info_path = repo_root / RELEASE_INFO_RELATIVE
    claim_registry_path = repo_root / CLAIM_REGISTRY_RELATIVE
    text = release_info_path.read_text(encoding="utf-8")
    registry_text = claim_registry_path.read_text(encoding="utf-8")
    registry = json.loads(registry_text)
    # Keep the registry's existing escaping so a bump changes one line, not every
    # non-ASCII character in the file.
    ascii_escaped = json.dumps(registry, indent=2, ensure_ascii=True) + "\n" == registry_text

    current_release_id = extract_macro(text, RELEASE_ID_MACRO)
    current_release_date = extract_macro(text, RELEASE_DATE_MACRO)

    if args.retag_urls_only:
        return retag_urls_only(repo_root, args.from_release_id, current_release_id, args.dry_run)

    next_release_id = args.release_id or increment_release_id(current_release_id)
    next_release_date = args.date or today_release_date()

    updated_text = replace_macro(text, RELEASE_ID_MACRO, next_release_id)
    updated_text = replace_macro(updated_text, RELEASE_DATE_MACRO, next_release_date)

    url_rewrites = {}
    if not args.keep_release_urls:
        url_rewrites = plan_release_url_rewrites(repo_root, current_release_id, next_release_id)

    if args.dry_run:
        print(f"{release_info_path}")
        print(f"release_id: {current_release_id} -> {next_release_id}")
        print(f"released_at: {current_release_date} -> {next_release_date}")
        print(f"claim_registry.release_id: {registry.get('release_id')} -> {next_release_id}")
        for path, count in sorted(url_rewrites.items()):
            print(f"release URLs: {path.relative_to(repo_root)}: {count}")
        return 0

    for path, _count in url_rewrites.items():
        path.write_text(
            rewrite_release_urls(path.read_text(encoding="utf-8"), current_release_id, next_release_id),
            encoding="utf-8",
        )

    registry["release_id"] = next_release_id
    release_info_path.write_text(updated_text, encoding="utf-8")
    claim_registry_path.write_text(
        json.dumps(registry, indent=2, ensure_ascii=ascii_escaped) + "\n",
        encoding="utf-8",
    )
    print(f"Updated {release_info_path}")
    print(f"Updated {claim_registry_path}")
    print(f"release_id: {current_release_id} -> {next_release_id}")
    print(f"released_at: {current_release_date} -> {next_release_date}")
    if url_rewrites:
        total = sum(url_rewrites.values())
        print(f"Rewrote {total} release-tagged repository URLs in {len(url_rewrites)} TeX sources")
    print(NEXT_STEP)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bump the shared paper release ID and release date in paper/release_info.tex for all current papers.",
    )
    parser.add_argument(
        "--release-id",
        help="Set an explicit release ID instead of incrementing the numeric suffix.",
    )
    parser.add_argument(
        "--date",
        help="Set an explicit release date instead of using today's local date.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the proposed release change without writing the file.",
    )
    parser.add_argument(
        "--retag-urls-only",
        action="store_true",
        help=(
            "Do not bump. Point repository URLs pinned to --from-release-id at the current "
            "release tag and stop (catch-up for sources that still cite an earlier release)."
        ),
    )
    parser.add_argument(
        "--from-release-id",
        help="Release tag the URLs currently cite; used with --retag-urls-only.",
    )
    parser.add_argument(
        "--keep-release-urls",
        action="store_true",
        help=(
            "Leave repository URLs that cite the current release tag (blob/<id>/, tree/<id>/) "
            "unchanged instead of pointing them at the new release."
        ),
    )
    return parser.parse_args()


def extract_macro(text: str, macro_name: str) -> str:
    pattern = re.compile(r"\\newcommand\{\\%s\}\{([^}]*)\}" % re.escape(macro_name))
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return match.group(1).strip()


def replace_macro(text: str, macro_name: str, value: str) -> str:
    pattern = re.compile(r"(\\newcommand\{\\%s\}\{)([^}]*)(\})" % re.escape(macro_name))
    updated_text, count = pattern.subn(r"\g<1>%s\g<3>" % value, text, count=1)
    if count != 1:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return updated_text


def retag_urls_only(repo_root: Path, from_release_id: str | None, current_release_id: str, dry_run: bool) -> int:
    if not from_release_id:
        raise SystemExit("--retag-urls-only needs --from-release-id <tag the URLs currently cite>")
    if from_release_id == current_release_id:
        raise SystemExit(f"URLs already cite the current release {current_release_id}")
    rewrites = plan_release_url_rewrites(repo_root, from_release_id, current_release_id)
    for path, count in sorted(rewrites.items()):
        print(f"release URLs: {path.relative_to(repo_root)}: {count}")
        if not dry_run:
            path.write_text(
                rewrite_release_urls(path.read_text(encoding="utf-8"), from_release_id, current_release_id),
                encoding="utf-8",
            )
    total = sum(rewrites.values())
    verb = "Would rewrite" if dry_run else "Rewrote"
    print(f"{verb} {total} repository URLs from {from_release_id} to {current_release_id} in {len(rewrites)} TeX sources")
    return 0


def release_url_pattern(release_id: str) -> re.Pattern[str]:
    """Match repository URLs pinned to ``release_id`` (``blob/<id>/`` or ``tree/<id>/``)."""
    return re.compile(
        r"(%s(?:blob|tree)/)%s(?=/)" % (re.escape(REPOSITORY_URL_FRAGMENT), re.escape(release_id))
    )


def rewrite_release_urls(text: str, current_release_id: str, next_release_id: str) -> str:
    return release_url_pattern(current_release_id).sub(r"\g<1>%s" % next_release_id, text)


def plan_release_url_rewrites(
    repo_root: Path, current_release_id: str, next_release_id: str
) -> dict[Path, int]:
    """Count the release-tagged repository URLs each TeX source would rewrite."""
    pattern = release_url_pattern(current_release_id)
    counts: dict[Path, int] = {}
    for directory in TEX_SOURCE_DIRS:
        base = repo_root / directory
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*.tex")):
            if "archive" in path.relative_to(repo_root).parts:
                continue
            found = len(pattern.findall(path.read_text(encoding="utf-8")))
            if found:
                counts[path] = found
    return counts


def increment_release_id(release_id: str) -> str:
    match = RELEASE_NUMBER_PATTERN.match(release_id)
    if not match:
        raise SystemExit(
            f"could not increment release id {release_id!r}; pass --release-id to set it explicitly"
        )
    prefix = match.group("prefix")
    number = int(match.group("number"))
    return f"{prefix}{number + 1}"


def today_release_date() -> str:
    today = datetime.now().astimezone().date()
    return f"{today.strftime('%B')} {today.day}, {today.year}"


if __name__ == "__main__":
    raise SystemExit(main())
