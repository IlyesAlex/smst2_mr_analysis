#%% Optional: in Spyder, run this once in the console if 3D rendering is flaky
# %matplotlib qt


#%% Imports
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_erosion


#%% Config
SCRIPT_DIR = Path.cwd()

COORD_PATH = SCRIPT_DIR / "../derivatives/figures/SI_figures/TableS2.xlsx"
TEMPLATE_PATH = SCRIPT_DIR / "../data/fmri/MNI_152_1mm.nii.gz"
MASK_PATH = SCRIPT_DIR / "../data/fmri/MNI_152_HC50_B.nii.gz"

OUTPUT_PATH = SCRIPT_DIR / "../derivatives/figures/Figure2C.png"
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

MAX_MASK_POINTS = 12000
SIZE_NORMAL = 55
SIZE_CURRENT = 120

# shell appearance
MASK_POINT_SIZE = 3
MASK_ALPHA_3D = 0.65
MASK_ALPHA_2D = 0.75   # slightly more transparent than before

# shrink factor for bottom-left 2D panel inside its quadrant
BOTTOM_LEFT_SHRINK = 0.92

# top-left and top-right 3D views
VIEWS_3D = [
    (18, 128),  # top-left
    (18, 40),   # top-right
]

AXIS_LABELSIZE = 14
AXIS_TICKSIZE = 12
LEGEND_FONTSIZE = 15


#%% Helper functions
def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Coordinate table not found: {path}")

    if path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep="\t", dtype=str)

    expected = {"id", "handle", "year", "raw", "x", "y", "z"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    return df


def to_float(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.replace(",", ".", regex=False)
        .replace({"": np.nan, "nan": np.nan, "None": np.nan})
        .astype(float)
    )


def get_mask_boundary_world(mask_img, max_points=12000):
    mask = mask_img.get_fdata() > 0

    eroded = binary_erosion(mask, iterations=1)
    boundary = mask & (~eroded)

    ijk = np.argwhere(boundary)
    if ijk.shape[0] == 0:
        raise ValueError("Mask appears empty.")

    if ijk.shape[0] > max_points:
        step = int(np.ceil(ijk.shape[0] / max_points))
        ijk = ijk[::step]

    xyz = nib.affines.apply_affine(mask_img.affine, ijk)
    return xyz


def get_mask_projection_xy(mask_img):
    """
    Create a filled 2D XY projection of the 3D mask and return:
    - proj_xy: 2D boolean array
    - extent: imshow extent in world coordinates [xmin, xmax, ymin, ymax]
    """
    mask = mask_img.get_fdata() > 0
    if not np.any(mask):
        raise ValueError("Mask appears empty.")

    # Project along z-axis -> XY occupancy
    proj_xy = np.any(mask, axis=2)

    nx, ny, _ = mask.shape

    # World coordinates of opposite XY corners
    p00 = nib.affines.apply_affine(mask_img.affine, [0, 0, 0])
    p11 = nib.affines.apply_affine(mask_img.affine, [nx - 1, ny - 1, 0])

    xmin, xmax = sorted([p00[0], p11[0]])
    ymin, ymax = sorted([p00[1], p11[1]])

    extent = [xmin, xmax, ymin, ymax]
    return proj_xy, extent


def set_equal_aspect_3d(ax, x, y, z, pad=5):
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    z_min, z_max = np.min(z), np.max(z)

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range) / 2.0

    x_mid = (x_max + x_min) / 2.0
    y_mid = (y_max + y_min) / 2.0
    z_mid = (z_max + z_min) / 2.0

    ax.set_xlim(x_mid - max_range - pad, x_mid + max_range + pad)
    ax.set_ylim(y_mid - max_range - pad, y_mid + max_range + pad)
    ax.set_zlim(z_mid - max_range - pad, z_mid + max_range + pad)


def set_equal_aspect_2d(ax, x, y, pad=8):
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)

    x_range = x_max - x_min
    y_range = y_max - y_min
    max_range = max(x_range, y_range) / 2.0

    x_mid = (x_max + x_min) / 2.0
    y_mid = (y_max + y_min) / 2.0

    ax.set_xlim(x_mid - max_range - pad, x_mid + max_range + pad)
    ax.set_ylim(y_mid - max_range - pad, y_mid + max_range + pad)
    ax.set_aspect("equal", adjustable="box")


def get_point_style(row, handle_to_color):
    marker = "x" if row["is_talairach"] else "o"

    if row["is_current"]:
        color = "red"
        size = SIZE_CURRENT
        lw = 2.0 if marker == "x" else 0.8
        edge = "black" if marker == "o" else None
        zorder = 10
    else:
        color = handle_to_color.get(row["handle"], "black")
        size = SIZE_NORMAL
        lw = 1.6 if marker == "x" else 0.5
        edge = "black" if marker == "o" else None
        zorder = 6

    return marker, color, size, lw, edge, zorder


def style_3d_axis(ax, mask_xyz):
    set_equal_aspect_3d(ax, mask_xyz[:, 0], mask_xyz[:, 1], mask_xyz[:, 2], pad=3)

    ax.set_xlabel("X (mm)", fontsize=AXIS_LABELSIZE, labelpad=8)
    ax.set_ylabel("Y (mm)", fontsize=AXIS_LABELSIZE, labelpad=8)
    ax.set_zlabel("Z (mm)", fontsize=AXIS_LABELSIZE, labelpad=6)

    ax.tick_params(axis="both", which="major", labelsize=AXIS_TICKSIZE, pad=2)

    ax.grid(False)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.fill = False
        axis.pane.set_edgecolor((1, 1, 1, 0))

    ax.set_box_aspect((1, 1, 1), zoom=1.08)

def style_2d_axis(ax, mask_xyz):
    ax.set_xlim(-46, 42)
    ax.set_ylim(-64, 22)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("X (mm)", fontsize=AXIS_LABELSIZE)
    ax.set_ylabel("Y (mm)", fontsize=AXIS_LABELSIZE)
    ax.tick_params(axis="both", which="major", labelsize=AXIS_TICKSIZE)

    ax.grid(False)

    ax.set_xticks(np.arange(-40, 41, 10))
    ax.set_yticks(np.arange(-60, 21, 10))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def draw_panel_3d(ax, df, mask_xyz, handle_to_color, elev, azim):
    ax.set_facecolor("white")

    # hippocampal shell
    ax.scatter(
        mask_xyz[:, 0],
        mask_xyz[:, 1],
        mask_xyz[:, 2],
        s=MASK_POINT_SIZE,
        c="lightgray",
        alpha=MASK_ALPHA_3D,
        depthshade=False,
        linewidths=0
    )

    # coordinates
    for _, row in df.iterrows():
        marker, color, size, lw, edge, zorder = get_point_style(row, handle_to_color)

        scatter_kwargs = dict(
            c=[color],
            s=size,
            marker=marker,
            linewidths=lw,
            depthshade=False,
            zorder=zorder
        )
        if marker == "o":
            scatter_kwargs["edgecolors"] = edge

        ax.scatter(row["x"], row["y"], row["z"], **scatter_kwargs)

    style_3d_axis(ax, mask_xyz)
    ax.view_init(elev=elev, azim=azim)


def draw_panel_xy(ax, df, mask_proj_xy, mask_extent_xy, mask_xyz, handle_to_color):
    ax.set_facecolor("white")

    # filled hippocampal mask in XY
    masked_proj = np.ma.masked_where(~mask_proj_xy, mask_proj_xy)
    ax.imshow(
        masked_proj.T,
        origin="lower",
        extent=mask_extent_xy,
        cmap=ListedColormap(["lightgray"]),
        alpha=MASK_ALPHA_2D,
        interpolation="nearest",
        zorder=1
    )

    # coordinates in XY only
    for _, row in df.iterrows():
        marker, color, size, lw, edge, zorder = get_point_style(row, handle_to_color)

        scatter_kwargs = dict(
            c=[color],
            s=size,
            marker=marker,
            linewidths=lw,
            zorder=zorder
        )
        if marker == "o":
            scatter_kwargs["edgecolors"] = edge

        ax.scatter(row["x"], row["y"], **scatter_kwargs)

    style_2d_axis(ax, mask_xyz)


def shrink_axis_in_place(ax, shrink=0.75):
    """
    Shrink an existing axes around its center while staying in the same quadrant.
    shrink=0.75 means 75% of the original width and height.
    """
    pos = ax.get_position()
    new_w = pos.width * shrink
    new_h = pos.height * shrink
    new_x = pos.x0 + (pos.width - new_w) / 2
    new_y = pos.y0 + (pos.height - new_h) / 2
    ax.set_position([new_x, new_y, new_w, new_h])


#%% Load coordinate table
df = load_table(COORD_PATH).copy()

df["id"] = df["id"].fillna("").astype(str).str.strip()
df["handle"] = df["handle"].fillna("").astype(str).str.strip()
df["raw"] = df["raw"].fillna("").astype(str)

df["x"] = to_float(df["x"])
df["y"] = to_float(df["y"])
df["z"] = to_float(df["z"])

df["is_current"] = df["id"].str.lower().eq("current")
df["is_talairach"] = df["raw"].str.contains("talairach", case=False, na=False)

df = df.dropna(subset=["x", "y", "z"]).reset_index(drop=True)

print(df[["id", "handle", "x", "y", "z", "is_current", "is_talairach"]])

if df.empty:
    raise ValueError("No valid coordinate rows found after parsing x/y/z.")


#%% Load images
template_img = nib.load(str(TEMPLATE_PATH))
mask_img = nib.load(str(MASK_PATH))

print("Template shape:", template_img.shape)
print("Mask shape:", mask_img.shape)

if template_img.shape != mask_img.shape:
    print("Warning: template and mask shapes differ.")

if not np.allclose(template_img.affine, mask_img.affine):
    print("Warning: template and mask affines differ.")

mask_xyz = get_mask_boundary_world(mask_img, max_points=MAX_MASK_POINTS)
mask_proj_xy, mask_extent_xy = get_mask_projection_xy(mask_img)

print("Boundary points used:", mask_xyz.shape[0])


#%% Set up colors
def extract_year_from_handle(handle: str) -> int:
    """
    Assumes the year is in the last 4 characters of the handle.
    Example: 'smith_2018' -> 2018
    """
    try:
        return int(str(handle)[-4:])
    except ValueError:
        return 9999  # fallback: push malformed handles to the end


non_current_handles = [h for h in df.loc[~df["is_current"], "handle"].unique() if h != ""]

# sort by year from the last 4 characters
non_current_handles = sorted(non_current_handles, key=extract_year_from_handle, reverse=True)

cmap = plt.get_cmap("tab20", max(len(non_current_handles), 1))
handle_to_color = {h: cmap(i) for i, h in enumerate(non_current_handles)}


#%% Create legend handles (no titles)
study_handles = []

if df["is_current"].any():
    study_handles.append(
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label="current",
            markerfacecolor="red",
            markeredgecolor="black",
            markersize=10
        )
    )

for handle in non_current_handles:
    study_handles.append(
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=handle,
            markerfacecolor=handle_to_color[handle],
            markeredgecolor="black",
            markersize=8
        )
    )

marker_handles = [
    Line2D(
        [0], [0],
        marker="o",
        color="black",
        linestyle="None",
        label="MNI-coordinates",
        markersize=7
    ),
    Line2D(
        [0], [0],
        marker="x",
        color="black",
        linestyle="None",
        label="visual inspection from Talairach-space",
        markersize=8
    ),
]


#%% Plot: 2x2 layout
fig = plt.figure(figsize=(14, 12), facecolor="white")

gs = fig.add_gridspec(
    2, 2,
    left=0.04,
    right=0.98,
    bottom=0.04,
    top=0.97,
    wspace=0.08,
    hspace=0.12
)

ax_tl = fig.add_subplot(gs[0, 0], projection="3d")
ax_tr = fig.add_subplot(gs[0, 1], projection="3d")
ax_bl = fig.add_subplot(gs[1, 0])      # 2D XY panel
ax_leg = fig.add_subplot(gs[1, 1])     # legend panel

# top-left and top-right: 3D with axes
for ax, (elev, azim) in zip([ax_tl, ax_tr], VIEWS_3D):
    draw_panel_3d(ax, df, mask_xyz, handle_to_color, elev, azim)

# bottom-left: top-down 2D XY with filled mask
draw_panel_xy(ax_bl, df, mask_proj_xy, mask_extent_xy, mask_xyz, handle_to_color)

# shrink bottom-left panel to 75% of its quadrant
shrink_axis_in_place(ax_bl, shrink=BOTTOM_LEFT_SHRINK)

# bottom-right: legend-only panel
ax_leg.axis("off")

legend_fontsize = LEGEND_FONTSIZE

legend1 = ax_leg.legend(
    handles=study_handles,
    loc="upper left",
    frameon=True,
    ncol=1,
    borderaxespad=0.0,
    bbox_to_anchor=(0.02, 0.98),
    fontsize=legend_fontsize,
    handlelength=1.6,
    labelspacing=0.9,
    borderpad=0.9,
    markerscale=1.15
)
ax_leg.add_artist(legend1)

legend2 = ax_leg.legend(
    handles=marker_handles,
    loc="lower left",
    frameon=True,
    ncol=1,
    borderaxespad=0.0,
    bbox_to_anchor=(0.02, 0.02),
    fontsize=legend_fontsize,
    handlelength=1.6,
    labelspacing=1.1,
    borderpad=0.9,
    markerscale=1.15
)


#%% Save and show
fig.canvas.draw()

fig.savefig(
    OUTPUT_PATH,
    dpi=300,
    facecolor="white",
    edgecolor="none",
    transparent=False
)

plt.show()

print(f"Saved figure to: {OUTPUT_PATH}")