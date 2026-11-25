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


"""
from tkinter import Tk, filedialog
from pathlib import Path
import numpy as np
import tifffile
from scipy.ndimage import distance_transform_edt


# User selects file to process
root = Tk()
root.withdraw()

selected_file = filedialog.askopenfilename(
    title="Select your input file"
)

root.destroy()

if not selected_file:
    raise SystemExit("No file selected.")


# Include selected file in parameter list, along with conversion specifications
FOLDER = Path(selected_file)
SPARSE = True
FRACTION = 0.02  # fraction of eligible pixels per class
MAX_PER_CLASS = 1000
BOUNDARY_PX = 1
SEED = 123


# Conversion function
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
        if mask.ndim != 2 or mask.shape != tifffile.imread(raw).shape:
            raise ValueError(f'Expected matching 2D images: {path}')
        if not np.isin(mask, [0, 1, 255]).all():
            raise ValueError(f'Not a binary mask: {path}')
        foreground = mask > 0
        labels = np.zeros(mask.shape, dtype=np.uint8)
        for region, label in [(foreground, 1), (~foreground, 2)]:
            if not sparse:
                labels[region] = label
                continue
            eligible = region & (distance_transform_edt(region) > BOUNDARY_PX)
            indices = np.flatnonzero(eligible)
            if not indices.size and region.any():
                print(f'Warning: no interior pixels for class {label} in {path.name}')
                indices = np.flatnonzero(region)
            if indices.size:
                n = min(MAX_PER_CLASS, max(1, round(FRACTION * indices.size)))
                labels.flat[rng.choice(indices, n, replace=False)] = label
        dest = output / ('Labels_' + raw.name.removeprefix('I_'))
        tifffile.imwrite(dest, labels, photometric='minisblack')
        print(dest.name, 'class counts:', np.bincount(labels.ravel(), minlength=3))


if __name__ == '__main__':
    convert(FOLDER, SPARSE)
