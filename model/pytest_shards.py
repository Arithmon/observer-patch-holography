"""Partition the complete pytest collection by module, preserving every node.

Load explicitly with ``-p model.pytest_shards``. Without both shard options,
the ordinary full-suite command remains unchanged. All jobs collect the same
suite before selecting modules, so a collection error fails every job.
"""
from pathlib import Path

import pytest


def selected_modules(modules, index, count):
    """A disjoint, exhaustive partition of sorted repository-relative paths."""
    if (type(index) is not int or type(count) is not int
            or count < 1 or not 0 <= index < count):
        raise ValueError("require 0 <= shard index < positive shard count")
    ordered = sorted(set(modules))
    if count > len(ordered):
        raise ValueError("shard count exceeds collected test modules")
    return frozenset(ordered[index::count])


def pytest_addoption(parser):
    group = parser.getgroup("model regression partitions")
    group.addoption("--model-shard-index", type=int, default=None)
    group.addoption("--model-shard-count", type=int, default=None)


def pytest_configure(config):
    index = config.getoption("model_shard_index")
    count = config.getoption("model_shard_count")
    if index is None and count is None:
        return
    if (index is None or count is None or count < 1 or not 0 <= index < count):
        raise pytest.UsageError("supply both model shard options with 0 <= index < count")
    if config.getoption("continue_on_collection_errors", default=False):
        raise pytest.UsageError("model shards require collection errors to fail closed")


@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    index = config.getoption("model_shard_index")
    count = config.getoption("model_shard_count")
    if index is None:
        return
    root = Path(config.rootpath)
    keys = [Path(item.path).relative_to(root).as_posix() for item in items]
    try:
        selected = selected_modules(keys, index, count)
    except ValueError as error:
        raise pytest.UsageError(str(error)) from error
    kept, removed = [], []
    for key, item in zip(keys, items, strict=True):
        (kept if key in selected else removed).append(item)
    config.hook.pytest_deselected(items=removed)
    items[:] = kept
    reporter = config.pluginmanager.getplugin("terminalreporter")
    if reporter is not None:
        reporter.write_line(
            f"model shard {index}/{count}: {len(kept)}/{len(keys)} collected tests, "
            f"{len(selected)}/{len(set(keys))} complete modules")
