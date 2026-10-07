#!/usr/bin/env python3
"""Audit FZ-13 CPL subsets and a base-LCDM display on official DESI chains.

This is a retrospective comparison.  DESI DR2 was public before this audit
and before any future FZ-13 freeze, so the output can diagnose conditional
branches but cannot score a prediction or confirm OPH.

The input files are the collaboration-produced Cobaya chains documented at
https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html.
Their hashes below are copied from the collaboration's v1.0 SHA-256 manifest.
No OPH code generated or refit these posterior samples.

The base-LambdaCDM companion calculation evaluates ``Lambda*l_P^2`` on each
weighted ``(H0, omegal)`` sample.  It therefore uses the chain covariance
directly and does not assume a correlation coefficient.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass, field
from decimal import Decimal
from fractions import Fraction
from numbers import Integral, Real
import re
from pathlib import Path
from typing import Iterable

import mpmath
from scipy.special import ndtri_exp


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SCRIPT_PATH = Path(__file__).resolve()
SOURCE_COBAYA_ROOT = (
    "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya"
)
SOURCE_MANIFEST = (
    "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
    "dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum"
)
SOURCE_MANIFEST_SHA256 = (
    "df78872aa8b2d3473a9e8de78f498180efd7cbcbeb18211ce4787fac52067ee5"
)

# Central constants used only to form the dimensionless retrospective display
# Lambda*l_P^2 = 3 Omega_Lambda (H0/c)^2 l_P^2 sample by sample.  The exact SI
# speed of light and the stated central conversion values are recorded in the
# receipt; their tiny uncertainties are not folded into the cosmology posterior.
SPEED_OF_LIGHT_M_S = 299_792_458.0
MPC_IN_M = 3.0856775814913673e22
PLANCK_LENGTH_M = 1.616255e-35
PLANCK_LENGTH_STANDARD_UNCERTAINTY_M = 0.000018e-35
CODATA_2022_PLANCK_LENGTH_SOURCE = "https://physics.nist.gov/cgi-bin/cuu/Value?plkl"

_CMB = (
    "planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_"
    "planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing"
)


DATASETS = {
    "DESI_DR2_BAO+CMB": {
        "slug": "cmb",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_{_CMB}",
        "sha256": [
            "c228de7bbaec19ddb22eec25c3dd7c40ef218976c020a7b55dd1c78dc3a638c5",
            "01bb30f43b3207d8575cc16354159fceff2ad4deffacc87964b0aef7f8e8ee44",
            "db12623c8c69c03ad219b28bbc517416fc941cbb1e717054344f30b4e43a4adc",
            "339312d9b7d3027147c38433fb4335cb9a12585776c4dfffc592c7c43af9a5ff",
        ],
    },
    "DESI_DR2_BAO+CMB+PantheonPlus": {
        "slug": "pantheon",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_pantheonplus_{_CMB}",
        "sha256": [
            "db81d299d59051ae5e1e8f67952320e08a1644a4a1f8d19267705aee32798877",
            "442390266f6a3bfe88e9c7cd8e29d1e7cff47fc8582f3a3af80a6b40b9f83939",
            "1c92edf523df34784633f6852f7564f3d953203ae939d5d2e8cd4c3fd6f580ae",
            "b17dc02689c3a30dd34a96d91ac35180006f17d10b82331396ef55ce999594af",
        ],
    },
    "DESI_DR2_BAO+CMB+Union3": {
        "slug": "union3",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_union3_{_CMB}",
        "sha256": [
            "f95d091203f615e1654ea9f4fea332db08fa66c59a5fe188b4915e67e1d73b11",
            "23aac58326257d207b8ee3699d21ed908e46238e42041d8c2ccf12fd85c8b839",
            "9e0998c5685e7552b94fd45e985dd137604f34692bebbc6fb922a1588669ab6a",
            "7710ababc2ab0c7fb5f3f67c8b790409703bce18707a2f99a1cd3fea531b4c20",
        ],
    },
    "DESI_DR2_BAO+CMB+DESY5": {
        "slug": "desy5",
        "model": "base_w_wa",
        "directory": f"desi-bao-all_desy5sn_{_CMB}",
        "sha256": [
            "8c783ebf283a205b7f569ce36a6694a32646ced5e28fd3b4683742733f6165e0",
            "cd4f2ff3a66aa92aceecd47c8452520c2de09d887f231fc0df033cd68597a886",
            "f7f8dbf28ff23d371e0b987b938d321adb920f1a09e97f746b3c0bb2d6247a01",
            "4a3607867d34890832f431c4d61947e5572a72ff1a407141f4b8cb8d3d7ec769",
        ],
    },
}

BASE_LCDM_NAME = "DESI_DR2_BAO+CMB_base_LCDM"
BASE_LCDM_DATASET = {
    "slug": "lcdm",
    "model": "base",
    "directory": f"desi-bao-all_{_CMB}",
    "sha256": [
        "00f3766f7a7b6370d21323886cd72869087b2b1346a04d729c8f3bc9e65ef698",
        "33b154eebdf4e9dca3b8f02ed2680120879d35c10b32fef42261a490104e1dc1",
        "d4717e7e5a13de851c86f24c87213faccef2b5f8747900274ab509d9dfa40aa2",
        "c827cd767a4864ca28aa15c902bda32004e803050d4be330e25aefddd78b5c36",
    ],
}

DOWNLOAD_SPECS = {**DATASETS, BASE_LCDM_NAME: BASE_LCDM_DATASET}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def producer_metadata() -> dict[str, str]:
    """Bind a receipt to the exact postprocessor that emitted it."""

    return {
        "script_path": SCRIPT_PATH.relative_to(REPO_ROOT).as_posix(),
        "script_sha256": sha256(SCRIPT_PATH),
    }


# Exact accumulation has no variance floor and is independent of row/merge order.
# Only final square roots and distribution functions use numerical arithmetic.
_MP = mpmath.mp.clone()
_MP.dps = 90
_REPORT_RTOL = Fraction(1, 10**12)
_DECIMAL_TOKEN = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def _fraction(value: object) -> Fraction:
    """Preserve supplied real scalars before any homogeneous float coercion."""
    if isinstance(value, bool):
        raise ValueError("expected a finite real scalar, not Boolean")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value))
    if isinstance(value, (Real, Decimal)):
        try:
            numerator, denominator = value.as_integer_ratio()
            return Fraction(int(numerator), int(denominator))
        except (AttributeError, ValueError, OverflowError) as exc:
            raise ValueError("expected a finite real scalar") from exc
    raise ValueError("expected a finite real scalar")


def _mp(value: Fraction):
    return _MP.mpf(value.numerator) / value.denominator


def _report(value, label: str) -> float:
    """Refuse erased or unresolved nonzero binary64 output (relative 1e-12)."""
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{label}: outside resolved binary64 range") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label}: outside resolved binary64 range")
    if isinstance(value, Fraction):
        error = abs(Fraction(result) - value)
        resolved = error <= abs(value) * _REPORT_RTOL
    else:
        resolved = abs(_MP.mpf(result) - value) <= abs(value) * _mp(_REPORT_RTOL)
    if not resolved:
        raise ValueError(f"{label}: outside resolved binary64 range")
    return result


def _sqrt(value: Fraction, label: str) -> float:
    return _report(_MP.sqrt(_mp(value)), label)


# Interpret the printed central conversion constants as decimal rationals.
_LAMBDA_FACTOR = (3 * (Fraction(1000) / Fraction(str(MPC_IN_M)) /
                      Fraction(str(SPEED_OF_LIGHT_M_S)))**2 *
                  Fraction(str(PLANCK_LENGTH_M))**2)


def _lambda_lp2(h0: Fraction, omega_lambda: Fraction) -> Fraction:
    if h0 <= 0 or not 0 < omega_lambda < 1:
        raise ValueError("H0 must be positive and OmegaLambda must lie in (0,1)")
    return _LAMBDA_FACTOR * omega_lambda * h0**2


def lambda_lp2_from_base_lcdm_sample(h0_km_s_mpc: float, omega_lambda: float) -> float:
    """Evaluate the declared SI conversion without rounded intermediate squares."""
    return _report(_lambda_lp2(_fraction(h0_km_s_mpc), _fraction(omega_lambda)),
                   "Lambda_lP2")


@dataclass
class WeightedMoments:
    """Exact weighted population moments, shared by CPL and base-LCDM."""
    dimension: int
    raw_rows: int = field(default=0, init=False)
    weight: Fraction = field(default=Fraction(0), init=False)
    weight_sq: Fraction = field(default=Fraction(0), init=False)
    sums: list[Fraction] = field(init=False)
    products: list[list[Fraction]] = field(init=False)

    def __post_init__(self) -> None:
        if self.dimension not in (2, 3):
            raise ValueError("posterior moments require two or three coordinates")
        self.sums = [Fraction(0) for _ in range(self.dimension)]
        self.products = [[Fraction(0) for _ in range(self.dimension)]
                         for _ in range(self.dimension)]

    def add(self, weight, values) -> None:
        weight = _fraction(weight)
        values = tuple(_fraction(value) for value in values)
        if weight <= 0:
            raise ValueError("chain weights must be positive")
        if len(values) != self.dimension:
            raise ValueError("posterior coordinate dimension mismatch")
        self.raw_rows += 1
        self.weight += weight
        self.weight_sq += weight**2
        for i, x in enumerate(values):
            self.sums[i] += weight * x
            for j in range(i, self.dimension):
                self.products[i][j] += weight * x * values[j]

    def merge(self, other: "WeightedMoments") -> None:
        if self.dimension != other.dimension:
            raise ValueError("posterior coordinate dimension mismatch")
        self.raw_rows += other.raw_rows
        self.weight += other.weight
        self.weight_sq += other.weight_sq
        self.sums = [a+b for a,b in zip(self.sums, other.sums, strict=True)]
        self.products = [[a+b for a,b in zip(row, other_row, strict=True)]
                         for row, other_row in zip(self.products, other.products, strict=True)]

    def mean(self, i: int) -> Fraction:
        if not self.weight:
            raise ValueError("empty chain")
        return self.sums[i] / self.weight

    def covariance(self, i: int, j: int) -> Fraction:
        i, j = sorted((i, j))
        return self.products[i][j] / self.weight - self.mean(i) * self.mean(j)

    def correlation(self, i: int, j: int) -> float | None:
        product = self.covariance(i, i) * self.covariance(j, j)
        cov = self.covariance(i, j)
        if not product:
            return None
        magnitude = _sqrt(cov**2 / product, "correlation")
        return -magnitude if cov < 0 else magnitude

    def counts(self) -> dict[str, float | int]:
        if not self.weight:
            raise ValueError("empty chain")
        return {
            "raw_rows": self.raw_rows,
            "expanded_posterior_weight": _report(self.weight, "weight"),
            "weight_concentration_ess_not_autocorrelation_corrected": _report(
                self.weight**2 / self.weight_sq, "weight concentration ESS"),
        }


def weighted_quantiles(samples, probabilities) -> dict[str, float]:
    """Step-CDF quantiles with exact weights and thresholds, including endpoints."""
    if not samples:
        raise ValueError("cannot compute quantiles of an empty sample")
    probabilities = tuple(_fraction(p) for p in probabilities)
    if any(not 0 <= p <= 1 for p in probabilities):
        raise ValueError("quantile probabilities must lie in [0,1]")
    if tuple(sorted(probabilities)) != probabilities:
        raise ValueError("quantile probabilities must be sorted")
    keys = [f"{float(p):.3f}" for p in probabilities]
    if len(set(keys)) != len(keys):
        raise ValueError("quantile probabilities have duplicate output keys")
    ordered = sorted((_fraction(value), _fraction(weight)) for value, weight in samples)
    if any(weight <= 0 for _, weight in ordered):
        raise ValueError("quantile samples require positive weights")
    total = sum(weight for _, weight in ordered)
    result = {}
    cumulative = Fraction(0)
    index = 0
    for value, weight in ordered:
        cumulative += weight
        while index < len(probabilities) and cumulative >= probabilities[index] * total:
            result[keys[index]] = _report(value, "quantile")
            index += 1
    return result


@dataclass
class BaseLCDMAccumulator:
    moments: WeightedMoments = field(default_factory=lambda: WeightedMoments(3))
    lambda_lp2_samples: list[tuple[Fraction, Fraction]] = field(default_factory=list)

    def add(self, weight, h0, omega_lambda) -> None:
        weight, h0, omega_lambda = map(_fraction, (weight, h0, omega_lambda))
        value = _lambda_lp2(h0, omega_lambda)
        self.moments.add(weight, (h0, omega_lambda, value))
        self.lambda_lp2_samples.append((value, weight))

    def merge(self, other: "BaseLCDMAccumulator") -> None:
        self.moments.merge(other.moments)
        self.lambda_lp2_samples.extend(other.lambda_lp2_samples)

    def summary(self) -> dict[str, object]:
        m = self.moments
        result = m.counts()
        for i, name in enumerate(("H0_km_s_Mpc", "OmegaLambda", "Lambda_lP2")):
            result[name] = {"weighted_mean": _report(m.mean(i), name + " mean"),
                            "weighted_std": _sqrt(m.covariance(i, i), name + " std")}
        result["H0_OmegaLambda_weighted_correlation"] = m.correlation(0, 1)
        result["Lambda_lP2"].update({
            "fractional_std_about_weighted_mean": _sqrt(
                m.covariance(2, 2) / m.mean(2)**2, "Lambda_lP2 fractional std"),
            "weighted_step_cdf_quantiles": weighted_quantiles(
                self.lambda_lp2_samples, tuple(map(Fraction, (".025", ".16", ".5", ".84", ".975")))),
        })
        return result


@dataclass
class Accumulator:
    moments: WeightedMoments = field(default_factory=lambda: WeightedMoments(2))
    monotone_weight: Fraction = Fraction(0)
    monotone_rows: int = 0
    w0_gt_neg_one_weight: Fraction = Fraction(0)
    wa_nonneg_weight: Fraction = Fraction(0)

    def add(self, weight, w0, wa) -> None:
        weight, w0, wa = map(_fraction, (weight, w0, wa))
        self.moments.add(weight, (w0, wa))
        # CPL is affine on [1/3,1]. Check the exact two endpoint inequalities.
        if w0 >= -1 and 3 * (w0 + 1) + 2 * wa >= 0:
            self.monotone_weight += weight
            self.monotone_rows += 1
        if w0 > -1:
            self.w0_gt_neg_one_weight += weight
        if wa >= 0:
            self.wa_nonneg_weight += weight

    def merge(self, other: "Accumulator") -> None:
        self.moments.merge(other.moments)
        for name in ("monotone_weight", "monotone_rows", "w0_gt_neg_one_weight", "wa_nonneg_weight"):
            setattr(self, name, getattr(self, name) + getattr(other, name))

    def summary(self) -> dict[str, float | int | None]:
        m = self.moments
        result = m.counts()
        for i, name in enumerate(("w0", "wa")):
            result[name + "_mean"] = _report(m.mean(i), name + " mean")
            result[name + "_std"] = _sqrt(m.covariance(i, i), name + " std")
        result.update({
            "w0_wa_covariance": _report(m.covariance(0, 1), "covariance"),
            "w0_wa_correlation": m.correlation(0, 1),
            "raw_rows_in_monotone_subset": self.monotone_rows,
        })
        for name, weight in (
            ("w_ge_minus_one_for_0_le_z_le_2", self.monotone_weight),
            ("capacity_loss_somewhere_for_0_le_z_le_2", m.weight - self.monotone_weight),
            ("w0_gt_minus_one", self.w0_gt_neg_one_weight),
            ("wa_nonnegative", self.wa_nonneg_weight),
        ):
            result["posterior_mass_" + name] = _report(weight / m.weight, name)
        return result


def _read_chain(path: Path, columns: tuple[str, ...], out):
    header = None
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split()
                    if len(set(header)) != len(header):
                        raise ValueError(f"{path}:{line_number}: duplicate column")
                    for column in columns:
                        if column not in header:
                            raise ValueError(f"{path}: required column absent: {column}")
                continue
            if header is None:
                raise ValueError(f"{path}:{line_number}: missing header")
            values = line.split()
            if len(values) != len(header):
                raise ValueError(f"{path}:{line_number}: {len(values)} values for {len(header)} columns")
            row = dict(zip(header, values, strict=True))
            try:
                if any(not _DECIMAL_TOKEN.fullmatch(row[c]) for c in columns):
                    raise ValueError("expected finite decimal columns")
                out.add(*(Fraction(row[c]) for c in columns))
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
    if header is None or out.moments.raw_rows == 0:
        raise ValueError(f"{path}: empty chain")
    return out


def read_chain(path: Path) -> Accumulator:
    return _read_chain(path, ("weight", "w", "wa"), Accumulator())


def read_base_lcdm_chain(path: Path) -> BaseLCDMAccumulator:
    return _read_chain(path, ("weight", "H0", "omegal"), BaseLCDMAccumulator())


def source_directory(spec: dict[str, object]) -> str:
    return f"{SOURCE_COBAYA_ROOT}/{spec['model']}/{spec['directory']}/"


def gaussian_fixed_point_diagnostic(accumulator: Accumulator) -> dict[str, object]:
    """Label an optional moment diagnostic without losing valid chain evidence.

    Invert the exact covariance, never the rounded standard deviations in the
    display. A singular covariance has no two-dimensional Gaussian diagnostic.
    """
    m = accumulator.moments
    v0, va, cov = m.covariance(0, 0), m.covariance(1, 1), m.covariance(0, 1)
    det = v0 * va - cov**2
    result = {
        "classification": "Gaussian moment summary; not official delta-chi2 or evidence",
        "status": "unavailable_singular_covariance",
        "mahalanobis_squared": None,
        "chi2_2dof_survival": None,
        "log_chi2_2dof_survival": None,
        "two_sided_normal_sigma_equivalent": None,
    }
    if det == 0:
        return result
    if det < 0:
        raise ArithmeticError("negative exact covariance determinant")
    d0, da = m.mean(0) + 1, m.mean(1)
    q = (va*d0**2 - 2*cov*d0*da + v0*da**2) / det
    try:
        result["mahalanobis_squared"] = _report(q, "Mahalanobis squared")
        result["log_chi2_2dof_survival"] = _report(-q/2, "log survival")
        # Small q needs expm1: subtracting a survival close to 1 erases sigma.
        sigma = (_MP.sqrt(2) * _MP.erfinv(-_MP.expm1(-_mp(q)/2)) if q <= 1
                 else -ndtri_exp(_report(-_mp(q)/2 - _MP.log(2), "log half-survival")))
        result["two_sided_normal_sigma_equivalent"] = _report(sigma, "normal sigma")
    except ValueError as exc:
        result.update(status="unavailable_binary64_range", reason=str(exc))
        return result
    try:
        # Avoid a huge exponential if the positive tail is certainly unresolved.
        if q > 1500:
            raise ValueError("survival: outside resolved binary64 range")
        result["chi2_2dof_survival"] = _report(_MP.exp(-_mp(q)/2), "survival")
        result["status"] = "available"
    except ValueError as exc:
        result.update(status="available_log_tail_only", reason=str(exc))
    return result


def audit_dataset(data_dir: Path, spec: dict[str, object]) -> dict[str, object]:
    combined = Accumulator()
    chains: list[dict[str, object]] = []
    expected_hashes = list(spec["sha256"])
    for index, expected in enumerate(expected_hashes, start=1):
        path = data_dir / f"{spec['slug']}_chain.{index}.txt"
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
        accumulator = read_chain(path)
        summary = accumulator.summary()
        chains.append(
            {
                "chain": index,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": actual,
                **summary,
            }
        )
        combined.merge(accumulator)

    summary = combined.summary()
    masses = [
        float(chain["posterior_mass_w_ge_minus_one_for_0_le_z_le_2"])
        for chain in chains
    ]
    raw_tail_rows = int(summary["raw_rows_in_monotone_subset"])
    return {
        "source_directory": source_directory(spec),
        "chains": chains,
        "combined": summary,
        "chain_range_for_monotone_subset_mass": [min(masses), max(masses)],
        "rare_tail_resolution": {
            "classification": (
                "resolved_in_all_four_chains"
                if min(int(chain["raw_rows_in_monotone_subset"]) for chain in chains)
                > 0
                else "at_least_one_chain_has_no_raw_tail_row"
            ),
            "combined_raw_tail_rows": raw_tail_rows,
            "warning": (
                "The weighted chain fraction is an empirical posterior mass under the DESI "
                "model, priors, and likelihoods. It is not a frequentist exclusion, a direct "
                "capacity measurement, or an OPH prediction score."
            ),
        },
        "fixed_capacity_point_gaussian_diagnostic": gaussian_fixed_point_diagnostic(combined),
    }


def audit_base_lcdm_dataset(
    data_dir: Path, spec: dict[str, object] = BASE_LCDM_DATASET
) -> dict[str, object]:
    """Compute the nonlinear Lambda*l_P^2 posterior from official samples."""

    combined = BaseLCDMAccumulator()
    chains: list[dict[str, object]] = []
    for index, expected in enumerate(spec["sha256"], start=1):
        path = data_dir / f"{spec['slug']}_chain.{index}.txt"
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
        accumulator = read_base_lcdm_chain(path)
        chains.append(
            {
                "chain": index,
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": actual,
                **accumulator.summary(),
            }
        )
        combined.merge(accumulator)
    return {
        "source_directory": source_directory(spec),
        "model_scope": (
            "flat base-LambdaCDM posterior under the official DESI DR2 BAO+CMB "
            "model, priors, and likelihoods"
        ),
        "sample_level_formula": (
            "Lambda*l_P^2 = 3*omegal*(H0*1000/Mpc_in_m/c)^2*l_P^2"
        ),
        "constants": {
            "speed_of_light_m_s_exact": SPEED_OF_LIGHT_M_S,
            "Mpc_in_m": MPC_IN_M,
            "Planck_length_m_CODATA_2022_central": PLANCK_LENGTH_M,
            "Planck_length_standard_uncertainty_m": (
                PLANCK_LENGTH_STANDARD_UNCERTAINTY_M
            ),
            "Planck_length_source": CODATA_2022_PLANCK_LENGTH_SOURCE,
            "constants_uncertainty_propagated": False,
        },
        "chains": chains,
        "combined": combined.summary(),
        "classification": (
            "sample-level retrospective display posterior; not an OPH prediction, "
            "model evidence, or model-independent measurement"
        ),
    }


def build_receipt(
    data_dir: Path, selected: Iterable[str] | None = None
) -> dict[str, object]:
    names = list(DATASETS) if selected is None else list(selected)
    unknown = sorted(set(names) - set(DATASETS))
    if unknown:
        raise ValueError(f"unknown datasets: {unknown}")
    return {
        "schema": "oph.official_desi_dr2_fz13_retrospective.v3",
        "producer": producer_metadata(),
        "arithmetic": {
            "chain_columns": "exact decimal values as published; no binary64 pre-rounding",
            "moments_and_subset_decisions": "exact rational population moments and CPL endpoints",
            "constants": "printed central SI decimal values; uncertainty not propagated",
            "numerical_reporting": "90-digit roots; binary64 fields require relative rounding error <= 1e-12",
            "gaussian_diagnostic": "exact moment inverse; numerical normal quantile; explicit unavailable status",
        },
        "source": {
            "publisher": "DESI Collaboration / DESI Data",
            "documentation": (
                "https://data.desi.lbl.gov/public/papers/y3/"
                "bao-cosmo-params/README.html"
            ),
            "official_sha256_manifest": SOURCE_MANIFEST,
            "official_sha256_manifest_sha256": SOURCE_MANIFEST_SHA256,
            "release_class": "official collaboration posterior chains",
        },
        "epistemic_status": {
            "retrospective_seen_data": True,
            "frozen_prediction_score": False,
            "oph_confirmation": False,
            "direct_capacity_measurement": False,
            "conditional_model_test_only": True,
            "notes": (
                "The CPL inequalities test only the conditional FZ-13 map after its density, "
                "closed-sector, and capacity-history premises. Agreement with (-1,0) is shared "
                "with LambdaCDM and earns no OPH-specific support."
            ),
        },
        "subset_definition": {
            "cpl": "w(a)=w0+wa(1-a)",
            "redshift_range": "0 <= z <= 2, equivalently 1/3 <= a <= 1",
            "monotone_capacity_condition": (
                "w(a)>=-1 on the full interval; because CPL is affine, this is exactly "
                "w0>=-1 and w0+(2/3)wa>=-1"
            ),
        },
        "datasets": {name: audit_dataset(data_dir, DATASETS[name]) for name in names},
        "base_lcdm_capacity_display": audit_base_lcdm_dataset(data_dir),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dataset", action="append", choices=sorted(DATASETS))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = build_receipt(args.data_dir, args.dataset)
    rendered = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
