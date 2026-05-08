#!/usr/bin/env bash
set -euo pipefail

BG="../data/fmri/sub-template_T1w_template_brain.nii.gz"

CA12_L="../data/fmri/masks/sub-template_FULLCA12_L.nii.gz"
CA12_R="../data/fmri/masks/sub-template_FULLCA12_R.nii.gz"
DGCA3_L="../data/fmri/masks/sub-template_FULLDGCA3_L.nii.gz"
DGCA3_R="../data/fmri/masks/sub-template_FULLDGCA3_R.nii.gz"
SUB_L="../data/fmri/masks/sub-template_FULLSUB_L.nii.gz"
SUB_R="../data/fmri/masks/sub-template_FULLSUB_R.nii.gz"

OUT="../derivatives/figures/SI_figures/Figure1SI_new.png"

# 30 coronal slices total: 89..118 inclusive
# (89..119 would be 31 slices)
Y_START=89
Y_END=118

# crop padding around bilateral MTL
PAD_X=8
PAD_Z=6

# per-panel render size
PANEL_W=420
PANEL_H=320

# gap between panels
GAP=12

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

mkdir -p "$(dirname "$OUT")"

for f in \
  "$BG" \
  "$CA12_L" "$CA12_R" \
  "$DGCA3_L" "$DGCA3_R" \
  "$SUB_L" "$SUB_R"
do
  [[ -f "$f" ]] || { echo "[ERROR] Missing file: $f" >&2; exit 1; }
done

for cmd in fslmaths fslstats fslroi fslval fsleyes pngappend; do
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "[ERROR] $cmd not in PATH" >&2
    exit 1
  }
done

# bilateral masks
fslmaths "$CA12_L"  -add "$CA12_R"  -bin "$tmpdir/CA12_bi.nii.gz"
fslmaths "$DGCA3_L" -add "$DGCA3_R" -bin "$tmpdir/DGCA3_bi.nii.gz"
fslmaths "$SUB_L"   -add "$SUB_R"   -bin "$tmpdir/SUB_bi.nii.gz"

# union for x/z crop box
fslmaths "$tmpdir/CA12_bi.nii.gz" -add "$tmpdir/DGCA3_bi.nii.gz" -add "$tmpdir/SUB_bi.nii.gz" -bin "$tmpdir/union.nii.gz"

read -r vox _mm3 < <(fslstats "$tmpdir/union.nii.gz" -V)
(( vox > 0 )) || { echo "[ERROR] union mask is empty" >&2; exit 1; }

read -r xmin xsize ymin ysize zmin zsize < <(
  fslstats "$tmpdir/union.nii.gz" -w | awk '{print $1, $2, $3, $4, $5, $6}'
)

dim1=$(fslval "$tmpdir/union.nii.gz" dim1)
dim2=$(fslval "$tmpdir/union.nii.gz" dim2)
dim3=$(fslval "$tmpdir/union.nii.gz" dim3)

xmax=$(( xmin + xsize - 1 ))
zmax=$(( zmin + zsize - 1 ))

x0=$(( xmin - PAD_X ))
x1=$(( xmax + PAD_X ))
z0=$(( zmin - PAD_Z ))
z1=$(( zmax + PAD_Z ))

(( x0 < 0 )) && x0=0
(( z0 < 0 )) && z0=0
(( x1 >= dim1 )) && x1=$(( dim1 - 1 ))
(( z1 >= dim3 )) && z1=$(( dim3 - 1 ))

xsize_crop=$(( x1 - x0 + 1 ))
ysize_crop=$dim2
zsize_crop=$(( z1 - z0 + 1 ))

(( xsize_crop > 0 && zsize_crop > 0 )) || {
  echo "[ERROR] Invalid crop size" >&2
  exit 1
}

(( Y_START >= 0 && Y_START < dim2 )) || { echo "[ERROR] Y_START out of range: $Y_START" >&2; exit 1; }
(( Y_END   >= 0 && Y_END   < dim2 )) || { echo "[ERROR] Y_END out of range: $Y_END" >&2; exit 1; }

# crop x/z only, keep full y so voxel indices stay valid
fslroi "$BG"                     "$tmpdir/bg_crop.nii.gz"    "$x0" "$xsize_crop" 0 "$ysize_crop" "$z0" "$zsize_crop"
fslroi "$tmpdir/CA12_bi.nii.gz"  "$tmpdir/ca12_crop.nii.gz"  "$x0" "$xsize_crop" 0 "$ysize_crop" "$z0" "$zsize_crop"
fslroi "$tmpdir/DGCA3_bi.nii.gz" "$tmpdir/dgca3_crop.nii.gz" "$x0" "$xsize_crop" 0 "$ysize_crop" "$z0" "$zsize_crop"
fslroi "$tmpdir/SUB_bi.nii.gz"   "$tmpdir/sub_crop.nii.gz"   "$x0" "$xsize_crop" 0 "$ysize_crop" "$z0" "$zsize_crop"

read -r dr_lo dr_hi < <(fslstats "$tmpdir/bg_crop.nii.gz" -P 2 -P 98)

echo "[INFO] Rendering 30 slices: ${Y_START}..${Y_END}"
echo "[INFO] Crop x: $x0 .. $x1"
echo "[INFO] Crop z: $z0 .. $z1"

# render individual slices
# ascending order so top-left becomes what used to be bottom-right
slice_pngs=()
idx=0
for y in $(seq "$Y_START" "$Y_END"); do
  png="$tmpdir/slice_$(printf '%02d' "$idx")_y${y}.png"

  fsleyes render \
    --outfile "$png" \
    --scene lightbox \
    --size "$PANEL_W" "$PANEL_H" \
    --zaxis y \
    --asVoxels \
    --zrange "$y" "$y" \
    --sliceSpacing 1 \
    --sampleSlices start \
    --labelSpace none \
    "$tmpdir/bg_crop.nii.gz" \
      -dr "$dr_lo" "$dr_hi" \
    "$tmpdir/ca12_crop.nii.gz" \
      -ot mask -mc 1 0 0 -o -w 2 -t 0.5 1.5 \
    "$tmpdir/dgca3_crop.nii.gz" \
      -ot mask -mc 0 0 1 -o -w 2 -t 0.5 1.5 \
    "$tmpdir/sub_crop.nii.gz" \
      -ot mask -mc 0 1 0 -o -w 2 -t 0.5 1.5

  slice_pngs+=( "$png" )
  idx=$(( idx + 1 ))
done

(( ${#slice_pngs[@]} == 30 )) || {
  echo "[ERROR] Expected 30 slice PNGs, got ${#slice_pngs[@]}" >&2
  exit 1
}

# build 5 rows x 6 columns
rows=()
for r in 0 1 2 3 4; do
  i=$(( r * 6 ))
  rowpng="$tmpdir/row_${r}.png"

  pngappend \
    "${slice_pngs[$((i+0))]}" + "$GAP" \
    "${slice_pngs[$((i+1))]}" + "$GAP" \
    "${slice_pngs[$((i+2))]}" + "$GAP" \
    "${slice_pngs[$((i+3))]}" + "$GAP" \
    "${slice_pngs[$((i+4))]}" + "$GAP" \
    "${slice_pngs[$((i+5))]}" \
    "$rowpng"

  rows+=( "$rowpng" )
done

pngappend \
  "${rows[0]}" - "$GAP" \
  "${rows[1]}" - "$GAP" \
  "${rows[2]}" - "$GAP" \
  "${rows[3]}" - "$GAP" \
  "${rows[4]}" \
  "$OUT"

echo "Saved: $OUT"
