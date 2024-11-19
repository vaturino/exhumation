import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import os


def compute_tm_counts(particles):
    """
    Calculate the count and percentage of particles with 1, 2, or 3 valid "tm" values.
    """
    # Count the number of non-NaN values across the "tm" columns
    particles["tm_count"] = particles[["tm_kin", "tm_dyn", "tm_trans"]].notna().sum(axis=1)

    # Determine the counts of particles with 1, 2, or 3 valid "tm" values
    tm_counts = particles["tm_count"].value_counts(normalize=True) * 100  # Convert to percentages
    one_tm = tm_counts.get(1, 0)  # Percentage of particles with 1 valid "tm"
    two_tm = tm_counts.get(2, 0)  # Percentage of particles with 2 valid "tm"
    return one_tm, two_tm


def main():
    # Load json file
    with open("models.json") as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/comparison/plots/velocity"
    os.makedirs(plot_loc, exist_ok=True)

    tests = ["velocity"]
    names = ["velocity_names"]

    # Initialize storage for plotting
    model_names = []
    one_tm_list = []
    two_tm_list = []

    for test in tests:
        if test in models:
            model_names = models[names[tests.index(test)]]
            for ind_m, m in enumerate(models[test]):
                text_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/txt_files"
                particles = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\s+")

                # Compute percentages of 1, 2, and 3 valid "tm" values
                one_tm, two_tm = compute_tm_counts(particles)
                one_tm_list.append(one_tm)
                two_tm_list.append(two_tm)

    # Plotting
    fig, ax = plt.subplots(figsize=(5, 5))
    bar_width = 0.3
    bar_positions = np.arange(len(model_names))

    # Stacked bar plot
    ax.bar(bar_positions, one_tm_list, color="cornflowerblue", label="1", width=bar_width)
    ax.bar(bar_positions, two_tm_list, color="lightcoral", label="2", bottom=one_tm_list, width=bar_width)

    # Labeling and formatting
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(model_names, rotation=45, ha="right")
    ax.legend(title="Number of\nStagnation\nintervals", fontsize=11, title_fontsize=11, loc = [0.6, 0.8])
    # Legend keys font size = 11
    # No spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_yticklabels([])
    ax.set_yticks([])

    # add percentage numbers in the bars
    for i, v in enumerate(one_tm_list):
        ax.text(i, v / 2, f"{v:.1f}%", color="black", ha="center", va="center", fontsize=11)
    for i, v in enumerate(two_tm_list):
        if v != 0:
            ax.text(i, v / 2 + one_tm_list[i], f"{v:.1f}%", color="black", ha="center", va="center", fontsize=11)


    plt.tight_layout()
    plt.savefig(f"{plot_loc}/tm_counts_comparison.eps", dpi=500)


if __name__ == "__main__":
    main()
