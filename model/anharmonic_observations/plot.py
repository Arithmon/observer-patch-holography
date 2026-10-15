"""Render all retained positive-duration comparisons; requires Matplotlib.

No model is selected, fitted or rerun here. The plot reads a completed receipt;
the panel choices, controls and displayed observation subsets are fixed below.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .data import FOLDER


CONTROLS = (
    ('primary', 'Executed pulse sequence', '#087e8b', '-'),
    ('frozen_control', 'Frozen preparation', '#6b7280', '--'),
    ('harmonic_control', 'Harmonic control', '#c55e24', '-.'),
    ('coherent_control', 'Zero-dephasing control', '#7851a9', ':'),
)


def render(result, output):
    """Render every prespecified CF point/time and the four fixed model choices."""
    runs = {run['name']: run['comparison'] for run in result['runs']}
    missing = {name for name, *_ in CONTROLS} - set(runs)
    if missing:
        raise ValueError(f'Missing prespecified plotting runs: {sorted(missing)}')
    primary = runs['primary']
    expected_rows = {'characteristic': 200, 'position': 507}
    for family, expected in expected_rows.items():
        if len(primary[family]['rows']) != expected:
            raise ValueError(f'Expected all {expected} positive-duration {family} means')
    if len(primary['derived_position']) != 39:
        raise ValueError('Expected all 39 positive-duration position diagnostics')

    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 12,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'oph-anharmonic-20261015'})
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.0))

    for axis, time in zip(axes[0], (2., 4.)):
        for component, color, marker, label in (
                ('real', '#087e8b', 'o', 'Real'),
                ('imag', '#7851a9', 's', 'Imaginary')):
            rows = sorted((row for row in primary['characteristic']['rows']
                           if row['time_ms'] == time and row['component'] == component),
                          key=lambda row: row['beta_imag'])
            if len(rows) != 50:
                raise ValueError(f'Expected all 50 {component} CF probes at {time:g} ms')
            probes = [row['beta_imag'] for row in rows]
            axis.plot(probes, [row['predicted_mean'] for row in rows], color=color,
                      linewidth=1.8, label=f'{label}: executed')
            axis.scatter(probes, [row['observed_mean'] for row in rows], color=color,
                         marker=marker, s=19, alpha=.75, linewidths=.5,
                         facecolors='none', label=f'{label}: measured', zorder=3)
        axis.axhline(0, color='#d1d5db', linewidth=.7, zorder=0)
        axis.set(title=f'Full characteristic-function scan at {time:g} ms',
                 xlabel=r'Imaginary displacement, Im $\beta$',
                 ylabel='Detector spin mean', xlim=(-.08, 5.08), ylim=(-1.08, 1.08))
        axis.grid(alpha=.15)
        axis.legend(loc='upper right', fontsize=8.5, frameon=True, framealpha=.92,
                    edgecolor='none', ncol=2, columnspacing=1.1, handlelength=1.8)

    position_axis = axes[1, 0]
    observed = sorted(primary['derived_position'], key=lambda row: row['time_ms'])
    position_axis.scatter([row['time_ms'] for row in observed],
                          [row['observed_position'] for row in observed],
                          s=23, color='#1f2937', label='Measured CF slope', zorder=4)
    for name, label, color, style in CONTROLS:
        rows = sorted(runs[name]['derived_position'], key=lambda row: row['time_ms'])
        if len(rows) != 39:
            raise ValueError(f'Expected all 39 position diagnostics for {name}')
        position_axis.plot([row['time_ms'] for row in rows],
                           [row['predicted_position'] for row in rows],
                           color=color, linestyle=style,
                           linewidth=2 if name == 'primary' else 1.5, label=label)
    position_axis.axhline(0, color='#d1d5db', linewidth=.7, zorder=0)
    position_axis.set(title='Position diagnostic from all 39 time slices',
                      xlabel='Published simulated duration (ms)',
                      ylabel='Position from CF slope (dimensionless)', xlim=(0, 16))
    position_axis.grid(alpha=.15)
    position_axis.legend(loc='upper center', bbox_to_anchor=(.5, -.22),
                         fontsize=8.5, ncol=2, frameon=False, columnspacing=1.3)

    error_axis = axes[1, 1]
    families = ('characteristic', 'position', 'combined')
    family_positions = np.arange(len(families))
    width = .19
    for index, (name, label, color, _) in enumerate(CONTROLS):
        values = [runs[name][family]['spin_mean_rmse'] for family in families]
        if not np.isfinite(values).all() or min(values) < 0:
            raise ValueError(f'Invalid descriptive RMSE for {name}')
        error_axis.bar(family_positions + width * (index - 1.5), values, width=width,
                       color=color, alpha=.85, label=label)
    error_axis.set_xticks(family_positions,
                         ['Full CF scans\n200 means', 'Time series\n507 means',
                          'Combined\n707 means'])
    error_axis.set(title='Descriptive error against every raw mean',
                   ylabel='Spin-mean RMSE', ylim=(0, None))
    error_axis.grid(axis='y', alpha=.15)
    error_axis.set_axisbelow(True)
    error_axis.legend(loc='upper center', bbox_to_anchor=(.5, -.22),
                      fontsize=8.5, ncol=2, frameon=False, columnspacing=1.3)

    fig.suptitle('Anharmonic evolution: executed observer readouts vs published measurements',
                 fontsize=16, y=.965)
    fig.text(.5, .925,
             'Preparation fitted only at zero duration; subsequent apparatus parameters fixed before scoring.',
             ha='center', fontsize=10, color='#374151')
    fig.text(.5, .035,
             'Retrospective comparison. Published apparatus, dephasing and clock are supplied.\n'
             'Aggregate means provide no authenticated error bars; position slopes reuse the same CF observations.',
             ha='center', fontsize=9, color='#374151', linespacing=1.6)
    fig.subplots_adjust(left=.075, right=.98, bottom=.205, top=.87,
                        wspace=.23, hspace=.39)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, metadata={'Software': 'Matplotlib'})
    plt.close(fig)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=FOLDER / 'results.json')
    parser.add_argument('--output', type=Path, default=FOLDER / 'comparison.png')
    args = parser.parse_args()
    render(json.loads(args.receipt.read_text()), args.output)


if __name__ == '__main__':
    main()
