"""
Compare hydrogen_excitation_rate arrays from two NetCDF files.
Usage: python compare_excitation_rates.py <old_file.nc> <new_file.nc>
"""

import sys
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(path):
    ds = nc.Dataset(path, "r")
    logT = ds.variables["logT"][:]
    rate = ds.variables["hydrogen_excitation_rate"][:]  # shape: (samples, start_level, final_level)
    ds.close()
    return logT, rate

def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_excitation_rates.py <old_file.nc> <new_file.nc>")
        sys.exit(1)

    old_path, new_path = sys.argv[1], sys.argv[2]
    logT_old, rate_old = load(old_path)
    logT_new, rate_new = load(new_path)

    if not np.allclose(logT_old, logT_new):
        print("WARNING: logT grids differ between files — comparison may be misleading.")

    # Absolute and relative differences
    abs_diff = np.abs(rate_new - rate_old)
    # Avoid divide-by-zero; mask where both are zero
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = np.where(
            np.abs(rate_old) > 0,
            abs_diff / np.abs(rate_old),
            np.nan
        )

    n_start = rate_old.shape[1]
    n_final = rate_old.shape[2]

    print(f"\n{'Pair':<12}  {'Max |rel diff|':>16}  {'Max |abs diff|':>16}  {'Mean |rel diff|':>16}")
    print("-" * 66)
    for i in range(1, n_start):
        for j in range(i + 1, n_final):
            rd = rel_diff[:, i, j]
            ad = abs_diff[:, i, j]
            valid = rd[~np.isnan(rd)]
            print(f"({i} -> {j})      "
                  f"  {np.nanmax(rd):16.4e}"
                  f"  {np.max(ad):16.4e}"
                  f"  {np.nanmean(rd):16.4e}")

    # --- Plot: one panel per (start, final) pair ---
    pairs = [(i, j) for i in range(1, n_start) for j in range(i + 1, n_final)
             if np.any(rate_old[:, i, j] != 0) or np.any(rate_new[:, i, j] != 0)]

    if not pairs:
        print("\nNo non-zero pairs found to plot.")
        return

    ncols = 2
    nrows = (len(pairs) + 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for idx, (i, j) in enumerate(pairs):
        ax = axes_flat[idx]
        ax.plot(logT_old, rate_old[:, i, j], label="old (hand-rolled)", lw=1.5)
        ax.plot(logT_new, rate_new[:, i, j], label="new (Boost)", lw=1.5, linestyle="--")
        ax.set_title(f"Level {i} → {j}")
        ax.set_xlabel("log₁₀ T")
        ax.set_ylabel("Excitation rate")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(len(pairs), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle("hydrogen_excitation_rate: old vs new", fontsize=14)
    plt.tight_layout()
    fig.savefig("excitation_rates_comparison.png", dpi=150)
    print("Saved excitation_rates_comparison.png")

    # --- Plot relative difference ---
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), squeeze=False)
    axes2_flat = axes2.flatten()

    for idx, (i, j) in enumerate(pairs):
        ax = axes2_flat[idx]
        ax.semilogy(logT_old, np.abs(rel_diff[:, i, j]), color="tab:red", lw=1.5)
        ax.set_title(f"Level {i} → {j}  |relative diff|")
        ax.set_xlabel("log₁₀ T")
        ax.set_ylabel("|Δ| / |old|")
        ax.grid(True, alpha=0.3)

    for idx in range(len(pairs), len(axes2_flat)):
        axes2_flat[idx].set_visible(False)

    fig2.suptitle("Relative difference (new vs old)", fontsize=14)
    plt.tight_layout()
    fig2.savefig("excitation_rates_reldiff.png", dpi=150)
    print("Saved excitation_rates_reldiff.png")
    plt.close("all")

if __name__ == "__main__":
    main()
