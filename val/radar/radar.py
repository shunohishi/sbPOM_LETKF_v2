import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Read CSV
# =========================
rmsd = pd.read_csv("RMSD.csv")
stat = pd.read_csv("Statistical_test.csv")

labels = rmsd["Observation"].values
N = len(labels)

# =========================
# Angle (CSV順・時計回り)
# =========================
angles = np.linspace(0, 2*np.pi, N, endpoint=False)
angles_closed = np.concatenate([angles, [angles[0]]])

# =========================
# Figure
# =========================
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# =========================
# Radius settings
# =========================
rmin, rmax = -50, 50
ax.set_ylim(rmin, rmax)
ax.set_yticks([-50, -25, 0, 25, 50])
ax.set_yticklabels(["-50", "-25", "0", "25", "50"], fontsize=11, color="0.35")
ax.set_rlabel_position(225)
ax.grid(False)

# =========================
# Background shading (+/-)
# =========================
theta_dense = np.linspace(0, 2*np.pi, 400)

# negative (blue)
ax.fill_between(
    theta_dense,
    rmin,
    0,
    color = "#f2f6fb",
    alpha=0.85,
    zorder=-10,
)

# positive (red)
ax.fill_between(
    theta_dense,
    0,
    rmax,
    color = "#fbf2ee",
    alpha=0.85,
    zorder=-10,
)

# =========================
# Circular grids
# =========================
theta_dense = np.linspace(0, 2*np.pi, 400)
for r in [-50, -25, 0, 25, 50]:
    lw = 2.8 if r == 0 else 1.0
    col = "black" if r == 0 else "0.7"
    ax.plot(theta_dense, np.full_like(theta_dense, r),
            color=col, lw=lw, zorder=0)

# =========================
# Radial guide lines
# =========================
for ang in angles:
    ax.plot([ang, ang], [rmin, rmax], color="0.75", lw=1, zorder=0)

# =========================
# Observation labels (robust outside placement)
# =========================
ax.set_xticks(angles)
ax.set_xticklabels([])

# base radius
label_r_base = rmax * 1.04

# labels that tend to intrude inward
extra_offset_labels = {
    "T (KEO)", "S (KEO)", "U (Papa)", "V (Papa)"
}

for ang, lab in zip(angles, labels):

    # extra offset only for problematic labels
    label_r = label_r_base * (1.10 if lab in extra_offset_labels else 1.00)

    # vertical alignment (force outward)
    if 0 <= ang < np.pi:
        va = "bottom"
    else:
        va = "top"

    # horizontal alignment
    if np.isclose(ang, 0) or np.isclose(ang, np.pi):
        ha = "center"
    elif 0 < ang < np.pi:
        ha = "left"
    else:
        ha = "right"

    ax.text(
        ang,
        label_r,
        lab,
        fontsize=16,
        ha=ha,
        va=va,
        clip_on=False,
        zorder=20,
    )

# =========================
# Models (LORAは描画しない)
# =========================
models = ["C-GLORS", "ORAS5", "GLORYS2V4"]
colors = {
    "C-GLORS": "#009E73",
    "ORAS5": "#E69F00",
    "GLORYS2V4": "#56B4E9",
}
linewidths = {
    "C-GLORS": 3.0,
    "ORAS5": 3.0,
    "GLORYS2V4": 3.0,
}
zorders = {
    "C-GLORS": 3,
    "ORAS5": 3,
    "GLORYS2V4": 3,
}

# Reference RMSD (LORA)
ref = rmsd["LORA-QG"].values

# =========================
# Plot RMSD ratio
# =========================
for m in models:
    ratio = (rmsd[m].values / ref - 1.0) * 100
    ratio_c = np.concatenate([ratio, [ratio[0]]])

    ax.plot(
        angles_closed, ratio_c,
        color=colors[m],
        lw=linewidths[m],
        label=m,
        zorder=zorders[m],
    )

    sig = stat[m].values
    for ang, val, s in zip(angles, ratio, sig):
        if s == 1:
            ax.scatter(
                ang, val, s=95,
                facecolors="white",
                edgecolors=colors[m],
                linewidths=2.3,
                zorder=10,
            )
        else:
            ax.scatter(
                ang, val, s=95,
                color=colors[m],
                zorder=10,
            )

# =========================
# Title & legend
# =========================
ax.set_title("RMSD ratio (%)", fontsize=24, pad=48)

handles, labels_leg = ax.get_legend_handles_labels()

# desired order
order = ["GLORYS2V4", "ORAS5", "C-GLORS"]

handles_ordered = [handles[labels_leg.index(o)] for o in order]
labels_ordered  = order

leg = ax.legend(
    handles_ordered,
    labels_ordered,
    loc="upper right",
    bbox_to_anchor=(1.30, 1.23),
    fontsize=15,
    frameon=True,
    facecolor="white",
    framealpha=0.9,
)

leg.get_frame().set_edgecolor("0.4")
leg.get_frame().set_linewidth(1.2)

# =========================
# Save & show
# =========================
plt.tight_layout()
plt.savefig("RMSD_ratio_radar.png", dpi=300, bbox_inches="tight")
plt.show()
