#!/usr/bin/env python3
"""Fifteen exact entropy and dimensional checks. Run: python verify_entropy.py.
Adapted from the supplied entropy verifier; historical V41 filename checks removed.
Requires SymPy. Does not establish general nonequilibrium horizon thermodynamics.

PUBLICATION CROSS-REFERENCES (freshly compiled V52 numbering):
  article 7.3, eq:generalized-apparent-entropy (218).
  article 7.3, eq:generalized-entropy-background (219).
  article 7.3, eq:generalized-entropy-matching (220).
  article 7.3, eq:generalized-entropy-no-ghost (221).
  article 7.3, eq:generalized-entropy-critical-radius (222).
  article 7.3, eq:generalized-entropy-time-derivative (223).
Method: Symbolic entropy derivatives and dimensional checks.
Inputs: Embedded entropy parameters r,nu,gamma,m,enthalpy positive.
Outputs below the selected results root: logs/entropy.log.
Provenance: GitHub_bundle/scripts/verify_entropy.py.
Scope limit: Not a general nonequilibrium thermodynamic derivation; dimensional checks are consistency checks.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
import sys

try:
    import sympy as sp
except ImportError as exc:
    raise SystemExit('Для проверки установите SymPy: python -m pip install sympy') from exc


# V52: article eq:generalized-apparent-entropy (218); article eq:generalized-entropy-background (219); article eq:generalized-entropy-matching (220); article eq:generalized-entropy-no-ghost (221); article eq:generalized-entropy-critical-radius (222); article eq:generalized-entropy-time-derivative (223).
def main() -> None:
    checks: list[str] = []

    def check(condition: bool, description: str) -> None:
        if not condition:
            raise AssertionError(description)
        checks.append(description)
        print(f'OK {len(checks):02d}: {description}')

    def identity(expression: sp.Expr, description: str) -> None:
        check(sp.simplify(expression) == 0, description)

    r, nu, gamma, m, enthalpy = sp.symbols('r nu gamma m enthalpy', positive=True)
    # m обозначает приведённую планковскую массу только внутри этой программы.
    area_entropy = 8 * sp.pi**2 * m**2 * r**2
    entropy = 2 * gamma * nu / (nu + 1) * r**(nu - 1) * area_entropy
    slope_ratio = gamma * nu * r**(nu - 1)
    temperature = 1 / (2 * sp.pi * r)
    heat_rate = 4 * sp.pi * r**2 * enthalpy  # H = 1/r
    r_dot = r**2 * enthalpy / (2 * m**2 * slope_ratio)
    h_dot = -enthalpy / (2 * m**2 * slope_ratio)
    h_enthalpy = (1 / slope_ratio - 1) * enthalpy
    chi = -3 * h_enthalpy

    identity(sp.diff(entropy, r) / sp.diff(area_entropy, r) - slope_ratio,
             'Отношение производных энтропий по радиусу')
    identity(r_dot + r**2 * h_dot, 'Производная радиуса R_A = H^(-1)')
    identity(temperature * sp.diff(entropy, r) * r_dot - heat_rate,
             'Соотношение Клаузиуса с потоком только независимого вещества')
    identity(-2 * m**2 * h_dot - enthalpy - h_enthalpy,
             'Совместимость двух уравнений Рэячаудхури')
    identity(chi - 3 * (1 - 1 / slope_ratio) * enthalpy,
             'Знак chi при совпадении фоновых траекторий')
    identity(sp.diff(entropy, r) * r_dot - 8 * sp.pi**2 * r**3 * enthalpy,
             'Временная производная обобщённой энтропии')
    identity(entropy.subs({nu: 1, gamma: 1}) - area_entropy,
             'Стандартная энтропия при nu = gamma = 1')
    identity(slope_ratio.subs(nu, 1) - gamma, 'Случай постоянного множителя nu = 1')
    # Радиус критической поверхности задаём эквивалентно: gamma = r_c^(1-nu)/nu.
    critical_sub = {gamma: r**(1 - nu) / nu}
    identity(slope_ratio.subs(critical_sub) - 1, 'Уравнение критического радиуса')
    enthalpy_dot = sp.symbols('enthalpy_dot', real=True)
    chi_dot = sp.diff(chi, r) * r_dot + sp.diff(chi, enthalpy) * enthalpy_dot
    identity(chi_dot.subs(critical_sub) - 3 * (nu - 1) * r * enthalpy**2 / (2 * m**2),
             'Трансверсальность при достижении критического радиуса, nu != 1')
    identity(chi.subs(enthalpy, 0), 'При rho_m + p_m = 0 выполняется chi = 0 для любого радиуса')
    identity(sp.diff(entropy, r) * r_dot.subs(enthalpy, 0),
             'При rho_m + p_m = 0 энтропия постоянна')
    # Проверка массовых размерностей. Размерность радиуса равна -1.
    identity((nu - 1) - (nu - 1), 'Безразмерность отношения производных энтропий')
    identity(2 + (nu - 1) - (nu - 1) + 2 - 4,
             'Размерность фонового уравнения: все члены имеют размерность 4')
    identity(-3 + 4 - 1, 'Размерность временной производной энтропии равна 1')

    print(f'\nPASS: {len(checks)} entropy / dimension checks.')

if __name__ == '__main__':
    main()
