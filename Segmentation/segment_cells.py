"""Run CP-SAM in the dedicated Cellpose environment."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np
import tifffile


def segment_cells(image_path=None, *, gpu=True, diameter=60, flow_threshold=2.0,
                  cellprob_threshold=-2, niter=10000, overwrite=False):
    from cellpose import models, io
    if image_path is None:
        from tkinter import Tk, filedialog
        root = Tk(); root.withdraw()
        try: image_path = filedialog.askopenfilename(title='Select I image')
        finally: root.destroy()
        if not image_path: return
    p = Path(image_path)
    if not p.name.startswith('I_') or p.suffix.lower() not in ('.tif','.tiff'):
        raise ValueError('Expected I_*.tif')
    b = p.stem[2:]
    labels = p.with_name(f'FPRM_Step2_CellLabels_{b}.tif')
    roi_base = p.with_name(f'FPRM_Step2_CellSegmentationOutlines_{b}')
    roi = Path(str(roi_base)+'_rois.zip')
    meta = p.with_name(f'FPRM_Step2_Metadata_{b}.json')
    for f in (labels,roi,meta):
        if f.exists() and not overwrite: raise FileExistsError(f)
    image = tifffile.imread(p)
    if image.ndim != 2: raise ValueError('Expected 2D image')
    model = models.CellposeModel(gpu=gpu)
    if gpu and not model.gpu: raise RuntimeError('GPU requested but unavailable')
    masks, _, _ = model.eval([image], diameter=diameter, flow_threshold=flow_threshold,
                             cellprob_threshold=cellprob_threshold, niter=niter)
    mask = np.asarray(masks[0])
    if mask.shape != image.shape or not np.issubdtype(mask.dtype,np.integer) or np.any(mask<0):
        raise ValueError('Invalid Cellpose labels')
    tifffile.imwrite(labels, mask.astype(np.uint32), photometric='minisblack')
    io.save_rois(mask, str(roi_base))
    if not roi.is_file(): raise RuntimeError('Expected ROI ZIP was not generated')
    meta.write_text(json.dumps(dict(schema_version='1.0', source_image=p.name,
        source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(), model='cpsam',
        cellpose_version=version('cellpose'), gpu=gpu, diameter=diameter,
        flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,niter=niter,
        processed_at=datetime.now(timezone.utc).isoformat(),
        outputs={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (labels,roi)}),indent=2))
    return labels,roi,meta

if __name__ == '__main__':
    segment_cells()
