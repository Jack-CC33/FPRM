"""Run CP-SAM in the dedicated Cellpose environment."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np
import tifffile


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


def segment_cells(image_path=None, *, metadata_path=None, model_path=None, gpu=True, diameter=60, flow_threshold=2.0,
                  cellprob_threshold=-2, niter=10000, overwrite=False):
    if image_path is None and metadata_path is None:
        metadata_path = select_experiment_metadata()
    initial_dir, experiment_reference = load_experiment_metadata(metadata_path)
    from cellpose import models, io
    if image_path is None:
        from tkinter import Tk, filedialog
        root = Tk(); root.withdraw()
        try: image_path = filedialog.askopenfilename(title='Select I image', initialdir=str(initial_dir) if initial_dir else '')
        finally: root.destroy()
        if not image_path: return
    p = Path(image_path)
    if not p.name.startswith('I_') or p.suffix.lower() not in ('.tif','.tiff'):
        raise ValueError('Expected I_*.tif')
    if initial_dir is not None and p.resolve().parent != initial_dir:
        raise ValueError('I image must be in the selected metadata directory.')
    b = p.stem[2:]
    labels = p.with_name(f'FPRM_Step2_CellLabels_{b}.tif')
    roi_base = p.with_name(f'FPRM_Step2_CellSegmentationOutlines_{b}')
    roi = Path(str(roi_base)+'_rois.zip')
    meta = p.with_name(f'FPRM_Step2_Metadata_{b}.json')
    for f in (labels,roi,meta):
        if f.exists() and not overwrite: raise FileExistsError(f)
    image = tifffile.imread(p)
    if image.ndim != 2: raise ValueError('Expected 2D image')
    model = models.CellposeModel(gpu=gpu, **({"pretrained_model": str(Path(model_path).resolve())} if model_path else {}))
    if gpu and not model.gpu: raise RuntimeError('GPU requested but unavailable')
    masks, _, _ = model.eval([image], diameter=diameter, flow_threshold=flow_threshold,
                             cellprob_threshold=cellprob_threshold, niter=niter)
    mask = np.asarray(masks[0])
    if mask.shape != image.shape or not np.issubdtype(mask.dtype,np.integer) or np.any(mask<0):
        raise ValueError('Invalid Cellpose labels')
    tifffile.imwrite(labels, mask.astype(np.uint32), photometric='minisblack')
    io.save_rois(mask, str(roi_base))
    if not roi.is_file(): raise RuntimeError('Expected ROI ZIP was not generated')
    meta.write_text(json.dumps(dict(schema_version='1.0', experiment_metadata=experiment_reference, source_image=p.name,
        source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(), model=Path(model_path).name if model_path else 'cpsam',
        model_path=str(Path(model_path).resolve()) if model_path else str(getattr(model, 'pretrained_model', 'cpsam')), 
        cellpose_version=version('cellpose'), gpu=gpu, diameter=diameter,
        flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,niter=niter,
        processed_at=datetime.now(timezone.utc).isoformat(),
        outputs={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (labels,roi)}),indent=2))
    return labels,roi,meta

if __name__ == '__main__':
    segment_cells()
