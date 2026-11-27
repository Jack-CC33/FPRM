"""Apply a trained two-class ilastik pixel classifier to a 2D TIFF."""

from pathlib import Path
from tempfile import TemporaryDirectory
import platform
import hashlib
import json
from datetime import datetime, timezone
import sys
import subprocess

import numpy as np
import tifffile


# USER SETTINGS
# Keep the local launcher path in ilastik_path.txt beside this script.
CONFIG_PATH = Path(__file__).resolve().with_name("ilastik_path.txt")
TIFF_SUFFIXES = (".tif", ".tiff")


def load_experiment_metadata(metadata_path):
    """Validate optional Step 1 metadata and return its portable directory and reference."""
    if metadata_path is None:
        return None, None
    import math
    path = Path(metadata_path).expanduser().resolve(strict=True)
    if path.name != "metadata.json":
        raise ValueError("Select the Step 1 metadata.json file.")
    raw = path.read_bytes()
    data = json.loads(raw.decode("utf-8-sig"))
    vit, free = data["SPvit"], data["SPfree"]
    if not all(type(v) in (int, float) and math.isfinite(v) for v in (vit, free)) or not 0 < vit < free:
        raise ValueError("Invalid FPRM calibration standards.")
    return path.parent, {"filename": path.name, "sha256": hashlib.sha256(raw).hexdigest()}


def select_experiment_metadata():
    """Offer metadata first; declining retains the original manual workflow."""
    from tkinter import Tk, filedialog, messagebox
    root = Tk()
    root.withdraw()
    try:
        if not messagebox.askyesno("FPRM setup", "Use Step 1 metadata.json to locate the experiment images?", parent=root):
            return None
        path = filedialog.askopenfilename(parent=root, title="Select experiment metadata.json", filetypes=[("JSON files", "*.json")])
        if not path:
            raise SystemExit("Metadata selection cancelled.")
        return path
    finally:
        root.destroy()


def is_wsl():
    """Detect Windows Subsystem for Linux, not ordinary Linux."""
    return (
        platform.system() == "Linux"
        and "microsoft" in platform.release().lower()
    )


def get_ilastik_path():
    """Read the local executable or launcher path (one line, no quotes)."""
    if not CONFIG_PATH.is_file():
        raise FileNotFoundError(
            f"Create {CONFIG_PATH} and enter your ilastik launcher path."
        )

    raw = CONFIG_PATH.read_text(encoding="utf-8-sig").strip()
    if not raw or len(raw.splitlines()) != 1:
        raise ValueError("ilastik_path.txt must contain exactly one path.")

    # WSL users may enter either a Windows drive path or a /mnt/c path.
    if is_wsl() and len(raw) >= 3 and raw[1] == ":" and raw[2] in "\\/":
        raw = subprocess.check_output(
            ["wslpath", "-u", raw], text=True
        ).strip()

    launcher = Path(raw).expanduser()
    if not launcher.is_absolute():
        launcher = CONFIG_PATH.parent / launcher

    if not launcher.is_file():
        raise FileNotFoundError(f"Ilastik launcher not found: {launcher}")

    return launcher


def ilastik_argument_path(path, launcher):
    """Translate file arguments only when WSL launches Windows ilastik."""
    path = Path(path).absolute()

    if is_wsl() and launcher.suffix.lower() == ".exe":
        return subprocess.check_output(
            ["wslpath", "-w", str(path)], text=True
        ).strip()

    return str(path)


def select_files(initial_dir=None):
    """Select the trained model and raw image using native file dialogs."""
    from tkinter import Tk, filedialog

    root = Tk()
    root.withdraw()
    try:
        model = filedialog.askopenfilename(
            parent=root,
            title="Select the Ilastik model file",
            filetypes=[("Ilastik models", "*.ilp")],
        )
        if not model:
            raise SystemExit("No model selected.")
        print(f"Selected model: {model}")

        image = filedialog.askopenfilename(
            parent=root,
            title="Select the raw FPRM image",
            initialdir=str(initial_dir) if initial_dir else "",
            filetypes=[("TIFF images", ("*.tif", "*.tiff", "*.TIF", "*.TIFF"))],
        )
        if not image:
            raise SystemExit("No image selected.")
        print(f"Selected image: {image}")
    finally:
        root.destroy()

    return Path(model), Path(image)


def segment(model_path=None, image_path=None, *, metadata_path=None, mito_label=1, output_path=None, overwrite=False):
    """Run ilastik and save a uint8 binary mitochondrial mask."""
    if model_path is None and image_path is None and metadata_path is None:
        metadata_path = select_experiment_metadata()
    initial_dir, experiment_reference = load_experiment_metadata(metadata_path)
    if model_path is None and image_path is None:
        model_path, image_path = select_files(initial_dir)
    elif model_path is None or image_path is None:
        raise ValueError("Supply both model_path and image_path, or neither.")

    model_path = Path(model_path).absolute()
    image_path = Path(image_path).absolute()
    if initial_dir is not None and image_path.resolve().parent != initial_dir:
        raise ValueError("I image must be in the selected metadata directory.")
    launcher = get_ilastik_path()

    for path in (model_path, image_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    if model_path.suffix.lower() != ".ilp":
        raise ValueError("Expected an ilastik .ilp project.")
    if not image_path.name.startswith("I_") or image_path.suffix.lower() not in TIFF_SUFFIXES:
        raise ValueError("Expected an I_*.tif or I_*.tiff input image.")
    if mito_label not in (1, 2):
        raise ValueError("mito_label must be 1 or 2.")

    image = tifffile.imread(image_path)
    if image.ndim != 2:
        raise ValueError("Expected a single-channel 2D image.")

    # Standardize final output to .tif, regardless of input extension.
    output = Path(output_path) if output_path else image_path.with_name("FPRM_Step3_MitoPrediction_" + image_path.stem[2:] + ".tif")
    metadata = output.with_name("FPRM_Step3_MitoPredictionMetadata_" + image_path.stem[2:] + ".json")
    if output.resolve() in (image_path.resolve(), model_path.resolve()):
        raise ValueError("Output cannot replace input.")
    if metadata.exists() and not overwrite:
        raise FileExistsError(metadata)
    if output.exists() and not overwrite:
        raise FileExistsError(output)

    # Keep temporary files beside the input so both WSL and Windows
    # can access them when the input is on a shared Windows drive.
    with TemporaryDirectory(prefix="mito_ilastik_", dir=image_path.parent) as temp:
        export_path = Path(temp) / "segmentation.tiff"

        command = [
            str(launcher),
            "--headless",
            f"--project={ilastik_argument_path(model_path, launcher)}",
            "--export_source=Simple Segmentation",
            "--output_format=tiff",
            f"--output_filename_format={ilastik_argument_path(export_path, launcher)}",
            ilastik_argument_path(image_path, launcher),
        ]

        print(f"Using ilastik: {launcher}")
        subprocess.run(command, check=True)

        # Ilastik versions may choose either .tif or .tiff.
        files = [
            p for p in Path(temp).iterdir()
            if p.is_file() and p.suffix.lower() in TIFF_SUFFIXES
        ]
        if len(files) != 1:
            raise RuntimeError(f"Expected one exported TIFF; found {files}")

        labels = tifffile.imread(files[0])
        # Remove an optional singleton channel axis, not spatial axes.
        if labels.ndim != 2:
            labels = labels.squeeze()
            if labels.ndim != 2 and labels.size == image.size:
                # Only restore singleton spatial dimensions.
                if tuple(s for s in image.shape if s != 1) == labels.shape:
                    labels = labels.reshape(image.shape)

        if labels.shape != image.shape or not np.isin(labels, [1, 2]).all():
            raise ValueError(
                "Unexpected export: check class labels, axes and export scaling."
            )

        mask = labels == mito_label

    tifffile.imwrite(output, mask.astype(np.uint8)*255, photometric="minisblack")
    metadata.write_text(json.dumps({
        "schema_version": "1.0", "experiment_metadata": experiment_reference, "source_image": image_path.name,
        "source_sha256": hashlib.sha256(image_path.read_bytes()).hexdigest(),
        "model": model_path.name, "model_path": str(model_path), "model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "mito_label": mito_label, "export_source": "Simple Segmentation",
        "prediction": output.name, "prediction_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "processed_at": datetime.now(timezone.utc).isoformat(), "ilastik_version": None,
        "python_version": sys.version, "numpy_version": np.__version__, "tifffile_version": tifffile.__version__
    }, indent=2))
    print(f"Saved {output} ({mask.sum():,} mitochondrial pixels)")
    return output


if __name__ == "__main__":
    segment()
