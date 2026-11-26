"""Generate binary compartment masks and a two-panel review."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import tifffile
from skimage.segmentation import find_boundaries


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


def generate_compartment_masks(image_path=None, cell_labels_path=None, mito_prediction_path=None, *, metadata_path=None, overwrite=False):
    if image_path is None and metadata_path is None:
        metadata_path = select_experiment_metadata()
    initial_dir, experiment_reference = load_experiment_metadata(metadata_path)
    if image_path is None:
        from tkinter import Tk, filedialog
        root=Tk(); root.withdraw()
        try: image_path=filedialog.askopenfilename(title='Select I image', initialdir=str(initial_dir) if initial_dir else '')
        finally: root.destroy()
        if not image_path: return
    p=Path(image_path)
    if not p.name.startswith('I_') or p.suffix.lower() not in ('.tif','.tiff'): raise ValueError('Expected I TIFF')
    if initial_dir is not None and p.resolve().parent != initial_dir:
        raise ValueError("I image must be in the selected metadata directory.")
    b=p.stem[2:]
    lp=Path(cell_labels_path or p.with_name(f'FPRM_Step2_CellLabels_{b}.tif'))
    mp=Path(mito_prediction_path or p.with_name(f'FPRM_Step3_MitoPrediction_{b}.tif'))
    image, labels, pred = map(tifffile.imread,(p,lp,mp))
    if image.ndim!=2 or labels.shape!=image.shape or pred.shape!=image.shape: raise ValueError('Shape mismatch')
    if not np.issubdtype(labels.dtype,np.integer) or np.any(labels<0): raise ValueError('Invalid instance labels')
    if not (np.isin(pred,[0,1]).all() or np.isin(pred,[0,255]).all()): raise ValueError('Invalid binary prediction')
    cell=labels>0; mito=(pred>0)&cell; cyto=cell&~mito
    if np.any(mito&cyto) or not np.array_equal(mito|cyto,cell): raise RuntimeError('Invalid partition')
    outputs={k:p.with_name(f'FPRM_Step3_{k}_{b}.tif') for k in ('CellMask','MitoMask','CytoMask','SegmentationReview')}
    meta=p.with_name(f'FPRM_Step3_Metadata_{b}.json')
    for f in (*outputs.values(),meta):
        if f.resolve() in [q.resolve() for q in (p,lp,mp)]: raise ValueError('Output overlaps input')
        if f.exists() and not overwrite: raise FileExistsError(f)
    upstream={}
    for name in (f'FPRM_Step2_Metadata_{b}.json',f'FPRM_Step3_MitoPredictionMetadata_{b}.json'):
        f=p.with_name(name)
        upstream[name]=json.loads(f.read_text()) if f.exists() else None
        if upstream[name] and upstream[name].get('source_sha256') != hashlib.sha256(p.read_bytes()).hexdigest(): raise ValueError('Upstream source hash mismatch')
    for k,m in zip(('CellMask','MitoMask','CytoMask'),(cell,mito,cyto)):
        tifffile.imwrite(outputs[k],m.astype(np.uint8)*255,photometric='minisblack')
    lo,hi=np.percentile(image,[0,99.5]); gray=np.rint(np.clip((image.astype(float)-lo)/(hi-lo if hi>lo else 1),0,1)*255).astype(np.uint8)
    rgb=np.repeat(gray[...,None],3,axis=2); overlay=rgb.copy()
    overlay[mito]=np.rint(.6*rgb[mito]+.4*np.array([0,255,255])).astype(np.uint8)
    overlay[find_boundaries(labels,mode='inner')]=[255,255,0]
    tifffile.imwrite(outputs['SegmentationReview'],np.concatenate([rgb,overlay],axis=1),photometric='rgb')
    meta.write_text(json.dumps(dict(schema_version='1.0',experiment_metadata=experiment_reference,basename=b,processed_at=datetime.now(timezone.utc).isoformat(),
        inputs={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (p,lp,mp)},
        outputs={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in outputs.values()},upstream=upstream,
        foreground_value=255,cytoplasm_definition='cell_minus_mito',manual_correction=False,
        review=dict(alpha=.4,display_limits=[float(lo),float(hi)],boundary_rgb=[255,255,0],mito_rgb=[0,255,255]),
        numpy_version=np.__version__,tifffile_version=tifffile.__version__),indent=2))
    return outputs

if __name__=='__main__': generate_compartment_masks()
