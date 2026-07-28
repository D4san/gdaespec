#!/usr/bin/env python3
"""Create consolidated plots for an isolated corrected G-DAE campaign."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "figure.dpi": 120,
    }
)

TESTS = [f"test_{i:02d}" for i in range(1, 6)]
CASES = [(0.00, 0.00), (0.01, 0.08), (0.08, 0.54), (0.26, 0.70)]
BRANCHES = ["phoenix", "sphinx"]
BRANCH_LABEL = {"phoenix": "PHOENIX", "sphinx": "SPHINX"}
BRANCH_COLOR = {"phoenix": "#2878B5", "sphinx": "#D95F02"}
TEST_COLORS = dict(zip(TESTS, ["#4477AA", "#CC6677", "#228833", "#EEAA33", "#AA3377"]))
PARAMS = ["log_CO2", "log_CH4", "log_O3", "log_H2O", "R_p_ref", "T"]
PARAM_LABELS = {
    "log_CO2": r"$\log_{10}(\mathrm{CO_2})$",
    "log_CH4": r"$\log_{10}(\mathrm{CH_4})$",
    "log_O3": r"$\log_{10}(\mathrm{O_3})$",
    "log_H2O": r"$\log_{10}(\mathrm{H_2O})$",
    "R_p_ref": r"$R_{p,\mathrm{ref}}\ (R_J)$",
    "T": r"$T\ (\mathrm{K})$",
}
EXPECTED = {
    "log_CO2": -3.0,
    "log_CH4": -8.0,
    "log_O3": -8.0,
    "log_H2O": -8.0,
    "R_p_ref": 0.0821,
    "T": 287.0,
}
NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PARAM_RE = re.compile(
    rf"^\s*(R_p_ref|T|log_H2O|log_CH4|log_CO2|log_O3)\s*=\s*"
    rf"({NUM})\s*\(\+({NUM})\)\s*\(-({NUM})\)"
)
CASE_RE = re.compile(r"(?P<spot>\d+\.\d+)spot-(?P<fac>\d+\.\d+)fac")


def case_text(spot: float, fac: float) -> str:
    return f"{spot:.2f}/{fac:.2f}"


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{stem}.png", out_dir / f"{stem}.pdf"]
    fig.savefig(paths[0], dpi=240, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def load_tables(base: Path) -> pd.DataFrame:
    metrics = pd.read_csv(base / "metrics.csv")
    times = pd.read_csv(base / "times.csv")
    keys = ["id", "branch", "f_spot", "f_fac", "strategy"]
    df = metrics.merge(times, on=keys, how="inner", validate="one_to_one")
    df["case"] = [case_text(s, f) for s, f in zip(df.f_spot, df.f_fac)]
    df["rmse_ppm"] = np.sqrt(df["MSE"]) * 1.0e6
    if len(df) != 35:
        raise ValueError(f"Expected 35 merged rows, found {len(df)}")
    return df


def plot_metrics(df: pd.DataFrame, out_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.8))
    fig.suptitle("Corrected G-DAE retrieval campaign — 35 cases", fontsize=18, y=0.995)

    ax = axes[0, 0]
    group_order = [
        (branch, spot, fac)
        for branch in BRANCHES
        for spot, fac in CASES
        if not (branch == "sphinx" and spot == 0.0 and fac == 0.0)
    ]
    for branch, spot, fac in group_order:
        sub = df[
            (df.branch == branch)
            & np.isclose(df.f_spot, spot)
            & np.isclose(df.f_fac, fac)
        ].sort_values("id")
        ax.plot(
            range(1, 6),
            sub.chi2_reduced,
            marker="o",
            ms=4,
            lw=1.25,
            alpha=0.78,
            color=BRANCH_COLOR[branch],
            label=BRANCH_LABEL[branch] if (spot, fac) == ((0.0, 0.0) if branch == "phoenix" else (0.01, 0.08)) else None,
        )
    ax.axhline(1.0, color="#444444", ls="--", lw=1.2, label=r"$\chi_r^2=1$")
    ax.set(xticks=range(1, 6), xlabel="Noise realization", ylabel=r"$\chi_r^2$", title="Reduced chi-square by realization")
    ax.set_ylim(0.25, 1.04)
    ax.legend(frameon=False, ncol=3, fontsize=9)

    ax = axes[0, 1]
    positions, labels = [], []
    pos = 0
    for branch in BRANCHES:
        for spot, fac in CASES:
            if branch == "sphinx" and spot == 0.0:
                continue
            sub = df[
                (df.branch == branch)
                & np.isclose(df.f_spot, spot)
                & np.isclose(df.f_fac, fac)
            ]
            ax.errorbar(
                pos,
                sub.chi2_reduced.mean(),
                yerr=sub.chi2_reduced.std(ddof=1),
                fmt="o",
                ms=7,
                capsize=4,
                color=BRANCH_COLOR[branch],
            )
            positions.append(pos)
            labels.append(f"{BRANCH_LABEL[branch][0]}\n{case_text(spot, fac)}")
            pos += 1
        pos += 0.55
    ax.axhline(1.0, color="#444444", ls="--", lw=1.2)
    ax.set(xticks=positions, xticklabels=labels, ylabel=r"$\chi_r^2$", title=r"Mean $\chi_r^2$ ± test scatter")
    ax.set_ylim(0.25, 1.04)

    ax = axes[1, 0]
    for branch, offset in [("phoenix", -0.13), ("sphinx", 0.13)]:
        for i, (spot, fac) in enumerate(CASES):
            if branch == "sphinx" and spot == 0.0:
                continue
            sub = df[
                (df.branch == branch)
                & np.isclose(df.f_spot, spot)
                & np.isclose(df.f_fac, fac)
            ]
            ax.errorbar(
                i + offset,
                sub.delta_time.mean(),
                yerr=sub.delta_time.std(ddof=1),
                fmt="o",
                ms=7,
                capsize=4,
                color=BRANCH_COLOR[branch],
                label=BRANCH_LABEL[branch] if i == (0 if branch == "phoenix" else 1) else None,
            )
    ax.set(
        xticks=range(4),
        xticklabels=[case_text(*c) for c in CASES],
        xlabel=r"$f_{\rm spot}/f_{\rm fac}$",
        ylabel="Time (min)",
        title="Retrieval runtime",
    )
    ax.legend(frameon=False)

    ax = axes[1, 1]
    for branch, offset in [("phoenix", -0.13), ("sphinx", 0.13)]:
        for i, (spot, fac) in enumerate(CASES):
            if branch == "sphinx" and spot == 0.0:
                continue
            sub = df[
                (df.branch == branch)
                & np.isclose(df.f_spot, spot)
                & np.isclose(df.f_fac, fac)
            ]
            ax.errorbar(
                i + offset,
                sub.rmse_ppm.mean(),
                yerr=sub.rmse_ppm.std(ddof=1),
                fmt="o",
                ms=7,
                capsize=4,
                color=BRANCH_COLOR[branch],
            )
    ax.set(
        xticks=range(4),
        xticklabels=[case_text(*c) for c in CASES],
        xlabel=r"$f_{\rm spot}/f_{\rm fac}$",
        ylabel="RMSE (ppm)",
        title="Spectral reconstruction error",
    )

    for ax in axes.ravel():
        ax.grid(True, ls=":", alpha=0.35)
        ax.tick_params(labelsize=9.5)
    fig.tight_layout()
    return save_figure(fig, out_dir, "gdae_corrected_campaign_metrics")


def parse_results(base: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(base.glob("test_*/**/retrievals/results/*_results.txt")):
        match = CASE_RE.search(path.name)
        if not match:
            continue
        branch = next((b for b in BRANCHES if b in path.parts), None)
        test_id = next((p for p in path.parts if p.startswith("test_")), None)
        if branch is None or test_id is None:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        block = text.split("1 σ constraints", 1)[-1].split("2 σ constraints", 1)[0]
        row = {
            "id": test_id,
            "branch": branch,
            "f_spot": float(match.group("spot")),
            "f_fac": float(match.group("fac")),
        }
        for line in block.splitlines():
            parsed = PARAM_RE.match(line)
            if parsed:
                name, value, plus, minus = parsed.groups()
                row[name] = float(value)
                row[f"{name}_plus"] = float(plus)
                row[f"{name}_minus"] = float(minus)
        if all(param in row for param in PARAMS):
            rows.append(row)
    result = pd.DataFrame(rows)
    if len(result) != 35:
        raise ValueError(f"Expected 35 parsed result files, found {len(result)}")
    return result


def plot_parameters(df: pd.DataFrame, out_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(3, 2, figsize=(13.2, 11.2), sharex=True)
    fig.suptitle("G-DAE atmospheric retrievals — posterior + realization uncertainty", fontsize=17, y=0.995)
    positions = {}
    labels = []
    pos = 0.0
    for branch in BRANCHES:
        for spot, fac in CASES:
            if branch == "sphinx" and spot == 0.0:
                continue
            positions[(branch, spot, fac)] = pos
            labels.append((pos, f"{BRANCH_LABEL[branch][0]}\n{case_text(spot, fac)}"))
            pos += 1.0
        pos += 0.6

    for ax, param in zip(axes.ravel(), PARAMS):
        for (branch, spot, fac), x in positions.items():
            sub = df[
                (df.branch == branch)
                & np.isclose(df.f_spot, spot)
                & np.isclose(df.f_fac, fac)
            ]
            center = float(sub[param].mean())
            scatter = float(sub[param].std(ddof=1))
            post_low = float(np.sqrt(np.mean(np.square(sub[f"{param}_minus"]))))
            post_up = float(np.sqrt(np.mean(np.square(sub[f"{param}_plus"]))))
            total_low = np.sqrt(scatter**2 + post_low**2)
            total_up = np.sqrt(scatter**2 + post_up**2)
            ax.errorbar(
                x,
                center,
                yerr=np.array([[total_low], [total_up]]),
                fmt="o",
                ms=6.5,
                capsize=3.5,
                color=BRANCH_COLOR[branch],
            )
        ax.axhline(EXPECTED[param], color="#333333", ls="--", lw=1.2, label="Injected value")
        ax.set_ylabel(PARAM_LABELS[param])
        ax.grid(True, ls=":", alpha=0.35)
        ax.set_xticks([x for x, _ in labels], [label for _, label in labels])
        ax.tick_params(axis="x", labelsize=8.5)
    axes[0, 0].legend(frameon=False, loc="best")
    axes[-1, 0].set_xlabel(r"Grid / $f_{\rm spot}/f_{\rm fac}$")
    axes[-1, 1].set_xlabel(r"Grid / $f_{\rm spot}/f_{\rm fac}$")
    fig.tight_layout()
    return save_figure(fig, out_dir, "gdae_corrected_retrieved_parameters")


def observation_path(base: Path, test_id: str, branch: str, spot: float, fac: float, recon: bool) -> Path:
    suffix = "_recon" if recon else ""
    branch_tag = "sphinx_" if branch == "sphinx" else ""
    return (
        base
        / test_id
        / branch
        / "observations"
        / f"pandexo_output_10transits_{branch_tag}fspot{spot:.2f}_ffac{fac:.2f}{suffix}.dat"
    )


def plot_reconstructions(base: Path, out_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.4), sharex=True)
    fig.suptitle("Corrected G-DAE spectra: noisy inputs and deterministic reconstructions", fontsize=17, y=0.995)
    for row, branch in enumerate(BRANCHES):
        for col, (spot, fac) in enumerate(CASES):
            ax = axes[row, col]
            if branch == "sphinx" and spot == 0.0:
                ax.axis("off")
                ax.text(0.5, 0.5, "No clean\nSPHINX case", ha="center", va="center", color="#666666")
                continue
            recon_stack, err_stack = [], []
            for test_id in TESTS:
                raw = np.loadtxt(observation_path(base, test_id, branch, spot, fac, False))
                rec = np.loadtxt(observation_path(base, test_id, branch, spot, fac, True))
                wl = rec[:, 0]
                ax.plot(raw[:, 0], raw[:, 2] * 1e6, color="#999999", lw=0.55, alpha=0.18)
                ax.plot(wl, rec[:, 2] * 1e6, color=TEST_COLORS[test_id], lw=0.95, alpha=0.76)
                recon_stack.append(rec[:, 2] * 1e6)
                err_stack.append(rec[:, 3] * 1e6)
            median = np.median(recon_stack, axis=0)
            median_err = np.median(err_stack, axis=0)
            ax.fill_between(wl, median - median_err, median + median_err, color=BRANCH_COLOR[branch], alpha=0.12)
            ax.plot(wl, median, color="#111111", lw=1.35, label="Median reconstruction")
            ax.set_title(f"{BRANCH_LABEL[branch]}  {case_text(spot, fac)}", fontsize=10.5)
            ax.grid(True, ls=":", alpha=0.28)
            if col == 0:
                ax.set_ylabel("Transit depth (ppm)")
            if row == 1:
                ax.set_xlabel(r"Wavelength ($\mu$m)")
    handles = [
        plt.Line2D([0], [0], color="#999999", lw=1.2, label="Noisy input"),
        plt.Line2D([0], [0], color="#111111", lw=1.5, label="Median reconstruction"),
        plt.Rectangle((0, 0), 1, 1, color=BRANCH_COLOR["phoenix"], alpha=0.16, label=r"Median reconstructed $1\sigma$"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return save_figure(fig, out_dir, "gdae_corrected_reconstruction_spectra")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    base = args.campaign_dir.resolve()
    out_dir = args.output_dir.resolve()

    outputs = []
    outputs.extend(plot_metrics(load_tables(base), out_dir))
    outputs.extend(plot_parameters(parse_results(base), out_dir))
    outputs.extend(plot_reconstructions(base, out_dir))
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
