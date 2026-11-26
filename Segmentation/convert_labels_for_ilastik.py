"""Convert paired binary TIFF masks to ilastik training labels.

Intention:
Ilastik is a very user-friendly software tool to help with semantic segmentation of diverse features within images, 
which uses Random Forest classifier to assign pixels to classes defined by a user.  For RAM usage consideration,
it is recommended to not over-label training images, but rather to adopt a 'sparse labeling' approach.  This script
enables a user with binary 0=background, 1=mitochondria (etc.) labeled TIFF images to convert these to sparsely-labeled
mask images with 0=unlabeled, 1=mitochondria, 2=non-mitochondria labeling.

Before this script:
I developed a python program that opens a graphical user interface that helps a user efficiently create binary segmentation
'ground truth' masks.  You can cycle through raw, representative images and create fully-labeled mask images with 0=background
and 1=foreground (in this case, mitochondria).

Run "roi_annotator.py," which may appear in this repository or may be separately provided.  Choose images to segment and draw
suitable masks.

This script:
Use the graphical user interface to select the FOLDER containing the raw images and corresponding binary masks.
The script will create a subfolder called "ilastik_labels" and write new TIFF images with the same names as the raw images,
but with the prefix "Labels_" and the "I_" prefix removed.  The new TIFF images will contain 0=unlabeled, 1=mitochondria, 
2=non-mitochondria values.

After this script:
Create a new pixel classification project in ilastik.
Under "Input Data", choose the raw images ("I_" prefixed) as input.
Under "Feature Selection", select desired (or all) parameters.
Under "Training", re-name the classes to "mitochondria" and "non-mitochondria."
      In the bottom-left corner, right-click on "labels" and choose "import," then find the "Labels_" prefixed file
        corresponding to the raw image you are currently training on.  In the dropdown just above "labels," select the next image and repeat.
Click "Live Update" to see the classifier's performance and make edits as needed.
When happy, save the .ilp model file to the FPRM>Segmentation>ilastik_models_mitochondria_segmentation folder associated with this repo.
"""

from tkinter import Tk, filedialog
from pathlib import Path
import numpy as np
import tifffile
from scipy.ndimage import distance_transform_edt


# User selects file to process
root = Tk()
root.withdraw()

folder_path = filedialog.askdirectory(title="Select the folder with masks")

if folder_path:
    print(f"Selected folder: {folder_path}")
else:
    raise SystemExit("No folder selected.")

root.destroy()


# Include selected file in parameter list, along with conversion specifications
FOLDER = Path(folder_path)
SPARSE = True
FRACTION = 0.03  # fraction of eligible pixels per class
MAX_PER_CLASS = 5000
SEED = 123

# Negative labels only: distance to the nearest mitochondrial pixel:
BOUNDARY_PX = 1
NEAR_PX = 10
FAR_PX = 30
BAND_WEIGHTS = (10, 3, 1)  # target counts: near, intermediate, far

# Negative labels only: intensity percentile bands, brightest first.
# Percentiles are calculated over all non-mitochondrial pixels in each image
INTENSITY_PERCENTILES = (50, 90)
INTENSITY_WEIGHTS = (10, 3, 1)  # bright, intermediate, dim


def sample_count(capacity):
    """Number of labels to sample from one class."""
    if capacity == 0:
        return 0
    return min(MAX_PER_CLASS, max(1, round(FRACTION * capacity)), capacity)

def allocate_counts(capacity, weights, n):
    """Allocate target counts, redistributing unavailable quotas."""
    capacity = np.asarray(capacity, dtype=int)
    weights = np.asarray(weights, dtype=float)
    if weights.shape != capacity.shape or not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError('Expected one positive finite weight per band')
    if not 0 <= n <= capacity.sum():
        raise ValueError('Requested count exceeds available pixels')
    counts = np.zeros_like(capacity)
    while counts.sum() < n:
        remaining = n - counts.sum()
        available = capacity - counts
        active_weights = weights * (available > 0)
        target = remaining * active_weights / active_weights.sum()
        extra = np.minimum(np.floor(target).astype(int), available)
        if not extra.any():
            extra[np.argmax(np.where(available > 0, target, -1))] = 1
        counts += extra
    return counts


def sample_mitochondria(region, rng):
    """Uniform random sampling: no intensity or distance preference."""
    indices = np.flatnonzero(region)
    return rng.choice(indices, size=sample_count(len(indices)), replace=False)


def sample_nonmitochondria(region, image, rng, report=False):
    """Allocate intensity quotas first, then distance quotas within each."""
    if not region.any():
        return np.empty(0, dtype=np.intp)
    if not (0 <= BOUNDARY_PX < NEAR_PX < FAR_PX):
        raise ValueError('Expected 0 <= BOUNDARY_PX < NEAR_PX < FAR_PX')
    low, high = INTENSITY_PERCENTILES
    if not 0 < low < high < 100:
        raise ValueError('Expected 0 < low < high < 100 percentiles')

    if region.all():
        distance = np.full(region.shape, np.inf)
    else:
        distance = distance_transform_edt(region)
    distance = distance.ravel()

    # Rank all negative pixels by raw intensity. Shuffle first so that
    # equal gray values at a percentile cutoff are split randomly.
    indices = rng.permutation(np.flatnonzero(region))
    values = image.ravel()[indices]
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite intensity in nonmitochondrial region')
    ranked = indices[np.argsort(values, kind='stable')]
    cut_low = round(len(ranked) * low / 100)
    cut_high = round(len(ranked) * high / 100)
    intensity_bands = [ranked[cut_high:], ranked[cut_low:cut_high], ranked[:cut_low]]

    # Retain the original negative boundary exclusion, if possible.
    eligible = distance[indices] > BOUNDARY_PX
    if eligible.any():
        intensity_bands = [b[distance[b] > BOUNDARY_PX] for b in intensity_bands]
    else:
        print('Warning: no negative interior pixels; including boundary pixels')

    capacity = [len(b) for b in intensity_bands]
    quotas = allocate_counts(capacity, INTENSITY_WEIGHTS, sample_count(sum(capacity)))
    samples = []
    joint_counts = []
    for band, quota in zip(intensity_bands, quotas):
        d = distance[band]
        distance_bands = [band[d <= NEAR_PX],
                          band[(d > NEAR_PX) & (d <= FAR_PX)],
                          band[d > FAR_PX]]
        counts = allocate_counts([len(b) for b in distance_bands], BAND_WEIGHTS, quota)
        joint_counts.append(counts)
        for candidates, count in zip(distance_bands, counts):
            samples.append(rng.choice(candidates, size=count, replace=False))
    if report:
        print('  Negative intensity counts (bright, middle, dim):', quotas.tolist())
        print('  Negative distance counts (near, middle, far):', np.sum(joint_counts, axis=0).tolist())
    return np.concatenate(samples)


def convert(folder, sparse=True):
    rng = np.random.default_rng(SEED)
    output = folder / 'ilastik_labels'
    output.mkdir(exist_ok=True)
    for path in sorted(folder.glob('I_*_mask.tif')):
        raw = path.with_name(path.stem.removesuffix('_mask') + '.tif')
        if not raw.is_file():
            print(f'Skipping unmatched mask: {path.name}')
            continue
        mask = tifffile.imread(path)
        image = tifffile.imread(raw)
        if mask.ndim != 2 or mask.shape != image.shape:
            raise ValueError(f'Expected matching 2D images: {path}')
        if not np.isin(mask, [0, 1, 255]).all():
            raise ValueError(f'Not a binary mask: {path}')
        foreground = mask > 0
        labels = np.zeros(mask.shape, dtype=np.uint8)
        if sparse:
            labels.flat[sample_mitochondria(foreground, rng)] = 1
            labels.flat[sample_nonmitochondria(~foreground, image, rng, report=True)] = 2
        else:
            labels[foreground] = 1
            labels[~foreground] = 2
        dest = output / ('Labels_' + raw.name.removeprefix('I_'))
        tifffile.imwrite(dest, labels, photometric='minisblack')
        print(dest.name, 'class counts:', np.bincount(labels.ravel(), minlength=3))

if __name__ == '__main__':
    convert(FOLDER, SPARSE)
