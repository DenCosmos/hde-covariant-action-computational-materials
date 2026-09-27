#!/usr/bin/env python3
"""Reconstruct the missing figure from eq:main-standard-cubic-scale.

No fitting or external data are used. On the standard family,
|chi| = 6 M_Pl^2 H^2 |c-1|/c. The plotted minimum is a bound from
individual cubic coefficients, not a physical unitarity threshold.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article A.11, fig:strong-coupling-standard (3).
  article A.11, eq:main-standard-cubic-scale (315).
  article A.11, eq:main-standard-cubic-numerical (317).
Method: Evaluate the two cubic operator estimates and their minimum.
Inputs: Embedded c grid and dimensionless normalization from original plotting script.
Outputs below the selected results root: figures/strong_coupling_scale_standard.csv, figures/strong_coupling_scale_standard.pdf.
Provenance: GitHub_bundle/scripts/plot_strong_coupling.py.
Scope limit: Operator estimates are not an independently determined physical unitarity threshold.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def normalized_bounds(c: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c = np.asarray(c, dtype=float)
    if np.any(c <= 0) or not np.all(np.isfinite(c)):
        raise ValueError('c must be finite and strictly positive')
    prefactor = (6.0 * np.abs(c-1.0) / c)**0.25
    clock = prefactor * 3.0**0.25 / 2.0
    c3 = prefactor * np.sqrt(9.0*np.sqrt(3.0)*np.abs(c-1.0)/(2.0*(2.0*c+7.0)))
    return clock, c3, np.minimum(clock,c3)

# V52: article fig:strong-coupling-standard (3); article eq:main-standard-cubic-scale (315); article eq:main-standard-cubic-numerical (317).
def main() -> None:
    base = Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=base/'results/current/figures')
    parser.add_argument('--data-output',type=Path,default=base/'results/current/figures')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    args.data_output.mkdir(parents=True,exist_ok=True)
    c=np.unique(np.r_[np.linspace(0.15,3.0,2201),0.55,1.0,1.2,25.0/16.0])
    clock,c3,bound=normalized_bounds(c)
    np.savetxt(args.data_output/'strong_coupling_scale_standard.csv',
       np.column_stack([c,clock,c3,bound]),delimiter=',',
       header='c,clock_over_sqrt_MH,C3_over_sqrt_MH,minimum_bound_over_sqrt_MH',comments='')
    fig,ax=plt.subplots(figsize=(7.2,4.35),layout='constrained')
    ax.plot(c,bound,linewidth=1.9,label='Minimum cubic-coefficient bound')
    point=float(normalized_bounds(np.array([1.2]))[2][0])
    ax.plot([1.2],[point],marker='o',linestyle='none')
    ax.annotate(r'$c=1.2$',xy=(1.2,point),xytext=(1.35,point-0.13),
                arrowprops=dict(arrowstyle='-'),fontsize=10)
    ax.axvline(1.0,linestyle=':',linewidth=0.8)
    ax.set_xlabel(r'$c$')
    ax.set_ylabel(r'Cubic-coefficient bound / $\sqrt{M_{\mathrm{Pl}}H}$')
    ax.set_xlim(0.15,3.0);ax.set_ylim(bottom=0)
    ax.tick_params(direction='in',top=True,right=True)
    fig.savefig(args.output/'strong_coupling_scale_standard.pdf')
    fig.savefig(args.data_output/'strong_coupling_scale_standard.png',dpi=160)
    plt.close(fig)
    near_one=6.0**0.25*np.sqrt(np.sqrt(3.0)/2.0)
    coefficient=float(np.sqrt(9*np.sqrt(3)*0.2/(2*(2*1.2+7))))
    report=(f'c=1.2 C3 dimensionless factor: {coefficient:.12g}\n'
            f'c=1.2 numerical eV bound for |chi|^1/4=1.9e-3 eV: {coefficient*1.9e-3:.12g}\n'
            f'near-c=1 coefficient: {near_one:.12g}\n'
            f'c=1.2 sound speed squared: {-(1.2+3)/(3*(1.2-1)):.12g}\n'
            'kinks: c=11/20 and c=25/16\n'
            'The zero at c=1 is the limiting coefficient bound, not a physical threshold.\n')
    (args.data_output/'strong_coupling_checks.txt').write_text(report)
    assert abs(coefficient-0.407)<0.0005
    assert abs(near_one-1.46)<0.005
    print(report)
if __name__=='__main__':main()
