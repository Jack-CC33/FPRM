"""Apply a trained two-class ilastik pixel classifier to a 2D TIFF."""
from tkinter import Tk, filedialog
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import tifffile
import subprocess

# User selects Ilastik model for segmentation
root = Tk()
root.withdraw()

model_path = filedialog.askopenfilename(title="Select the Ilastik model file", filetypes=[("Ilastik models", "*.ilp")])

if model_path:
    print(f"Selected model: {model_path}")
else:
    raise SystemExit("No model selected.")

root.destroy()

# User selects raw FPRM image to segment
root = Tk()
root.withdraw()

image_path = filedialog.askopenfilename(title="Select the raw FPRM image", filetypes=[("TIFF images", "*.tif")])

if image_path:
    print(f"Selected image: {image_path}")
else:
    raise SystemExit("No image selected.")

root.destroy()



# Change these paths for your computer.
MODEL_PATH = Path(model_path)
IMAGE_PATH = Path(image_path)
# None uses Easilastik's automatic executable discovery.
ILASTIK_PATH = Path(r"/mnt/c/Program Files/ilastik-1.4.2-gpu/ilastik.exe")
# Optional aligned Cellpose instance-label TIFF.
CELL_MASK_PATH = None
MITO_LABEL = 1
OVERWRITE = False

# Function needed to convert WSL paths to Windows paths for subprocess calls
# I run this python script in WSL, but my ilastik installation is in Windows.
# YOU MAY NEED TO MODIFY THIS SCRIPT IF USING A DIFFERENT SETUP.
def windows_path(path):
    """Convert a WSL filesystem path to a Windows path."""
    return subprocess.check_output(
        ["wslpath", "-w", str(path)],
        text=True
    ).strip()


def segment():
    from easilastik import run_ilastik

    for path in (MODEL_PATH, IMAGE_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not IMAGE_PATH.name.startswith('I_') or IMAGE_PATH.suffix.lower() != '.tif':
        raise ValueError('Expected an I_*.tif input image')
    image = tifffile.imread(IMAGE_PATH)
    if image.ndim != 2:
        raise ValueError('Expected a single-channel 2D image')
    output = IMAGE_PATH.with_name('MitoMask_' + IMAGE_PATH.name[2:])
    if output.exists() and not OVERWRITE:
        raise FileExistsError(output)

    # Export integer class labels, not a colored visualization.
    # A temporary directory prevents accidentally reading stale results.
    with TemporaryDirectory(prefix="mito_ilastik_", dir=IMAGE_PATH.parent) as temp:
        export_path = Path(temp) / "segmentation.tif"

        command = [
            str(ILASTIK_PATH),
            "--headless",
            f"--project={windows_path(MODEL_PATH)}",
            "--export_source=Simple Segmentation",
            "--output_format=tiff",
            f"--output_filename_format={windows_path(export_path)}",
            windows_path(IMAGE_PATH)
        ]

        subprocess.run(command, check=True)

        if not export_path.is_file():
            raise FileNotFoundError(
                f"Ilastik did not produce {export_path}"
            )

        labels = tifffile.imread(export_path).squeeze()

        if (
            labels.shape != image.shape
            or not np.isin(labels, [1, 2]).all()
        ):
            raise ValueError(
                "Unexpected export: check class labels, "
                "axes and export scaling"
            )

        mask = labels == MITO_LABEL

    if CELL_MASK_PATH is not None:
        cells = tifffile.imread(CELL_MASK_PATH)
        if cells.shape != image.shape:
            raise ValueError('Cellpose mask and image shapes differ')
        mask &= cells > 0
    tifffile.imwrite(output, mask.astype(np.uint8), photometric='minisblack')
    print(f'Saved {output} ({mask.sum():,} mitochondrial pixels)')
    return output


if __name__ == '__main__':
    segment()
