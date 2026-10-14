"""Render the retained comparison; requires optional Matplotlib."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .data import FOLDER


def main():
    result = json.loads((FOLDER/'results.json').read_text())
    comparison = result['heldout_comparison']
    rows = comparison['rows']
    plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'oph-readout-20261014'})
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5), gridspec_kw={'width_ratios': [1, 1.15]})
    for part, color, label in (('real', '#087e8b', 'Real component'),
                               ('imag', '#c55e24', 'Imaginary component')):
        axes[0].scatter([r['value']['prediction'][part] for r in rows],
                        [r['value']['observation'][part] for r in rows],
                        s=24, color=color, alpha=.7, label=label)
    axes[0].plot([-.5, .85], [-.5, .85], color='#4b5563', linestyle='--', linewidth=1)
    axes[0].set(xlabel='Prediction from the two measured axes', ylabel='Published spin mean',
                xlim=(-.5, .85), ylim=(-.5, .85), title='100 withheld complex responses')
    axes[0].legend(frameon=False, loc='upper left')
    axes[0].grid(alpha=.15)
    errors = np.empty((10, 10))
    for row in rows:
        i, j = row['index']
        errors[j-1, i-1] = row['value']['absolute_complex_error']
    image = axes[1].imshow(errors, origin='lower', cmap='magma', vmin=0, vmax=.25,
                            extent=(.5, 10.5, .5, 10.5))
    axes[1].set(xlabel='Real-displacement index', ylabel='Imaginary-displacement index',
                title='Absolute complex residual')
    fig.colorbar(image, ax=axes[1], fraction=.048, pad=.04)
    fig.suptitle('Published oscillator readout: axes calibrate, off-axis values test', fontsize=13)
    fig.text(.5, .02,
             'Per-component RMSE = 0.07639. Retrospective, postselected preparation/readout comparison; no dynamics tested.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .055, 1, .94))
    fig.savefig(FOLDER/'comparison.png', dpi=180, metadata={'Software': 'Matplotlib'})
    plt.close(fig)


if __name__ == '__main__':
    main()
