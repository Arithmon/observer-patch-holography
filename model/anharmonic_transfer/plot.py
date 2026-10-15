"""Plot all five fixed transfer tests and their descriptive errors."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FOLDER = Path(__file__).resolve().parent


def main():
    receipt = json.loads((FOLDER/'results.json').read_text())
    prediction = receipt['prediction']
    runs = {r['name']: r for r in prediction['runs']}
    values = [r['observed_position_estimator'] for r in receipt['scores']['primary']['rows']]
    values += [v['reported_zero_time_estimator'] for v in prediction['preparations'].values()]
    values += [r['predicted_position_estimator'] for run in runs.values() for r in run['predictions']]
    vertical_extent = max(1.8, 1.08*max(abs(v) for v in values))
    fig, axes = plt.subplots(3, 2, figsize=(13, 11), layout='constrained')
    descriptions = ['symmetric / left preparation', 'mild asymmetry / left preparation',
                    'mild asymmetry / right preparation', 'mild asymmetry / centred preparation',
                    'stronger asymmetry / left preparation']
    styles = [('primary', 'Transferred preparation', '#1565a8', '-'),
              ('preparation_adjusted', 'Zero-time preparation adjusted', '#d06b18', '--'),
              ('instantaneous_carrier_control', 'Instantaneous carrier control', '#707070', ':')]
    for ax, case, description in zip(axes.flat, 'abcde', descriptions, strict=False):
        measured = [r for r in receipt['scores']['primary']['rows'] if r['case'] == case]
        zero = prediction['preparations'][case]['reported_zero_time_estimator']
        ax.scatter([0]+[r['time_ms'] for r in measured], [zero]+[r['observed_position_estimator'] for r in measured],
                   s=21, color='#171717', zorder=5, label='Published derived measurements')
        for name, label, color, style in styles:
            rows = [r for r in runs[name]['zero_time_predictions']+runs[name]['predictions'] if r['case'] == case]
            ax.plot([r['time_ms'] for r in rows], [r['predicted_position_estimator'] for r in rows],
                    color=color, linestyle=style, linewidth=1.7, label=label)
        ax.set_title(f'({case}) {description}', loc='left', fontsize=11)
        ax.set(xlabel='Published simulated duration (ms)', ylabel='Two-point position estimator', xlim=(-.2, 10.2),
               ylim=(-vertical_extent, vertical_extent))
        ax.axhline(0, color='#bbbbbb', linewidth=.7)
        ax.grid(alpha=.14)
    bar = axes.flat[-1]
    names = ['primary', 'preparation_adjusted', 'instantaneous_carrier_control', 'harmonic_control', 'frozen_control']
    labels = ['Transferred preparation', 'Zero-time preparation adjusted', 'Instantaneous carrier', 'Harmonic control', 'Frozen control']
    errors = [receipt['scores'][name]['metrics']['combined']['position_estimator_rmse'] for name in names]
    bar.barh(np.arange(len(names)), errors, color=['#1565a8', '#d06b18', '#707070', '#999999', '#bbbbbb'])
    bar.set_yticks(np.arange(len(names)), labels, fontsize=9)
    bar.invert_yaxis()
    bar.set_xlabel('RMSE: all 125 positive-time derived observations')
    bar.set_title('Fixed comparisons; lower error is better', loc='left', fontsize=11)
    bar.set_xlim(0, max(errors)*1.22)
    for i, value in enumerate(errors):
        bar.text(value+max(errors)*.025, i, f'{value:.3f}', va='center', fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=2, fontsize=10, frameon=False)
    fig.suptitle('Changed-potential test: fixed execution versus published observations\n'
                 'Single-mode apparatus model · no positive-time calibration · descriptive comparison', fontsize=15)
    fig.savefig(FOLDER/'comparison.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
