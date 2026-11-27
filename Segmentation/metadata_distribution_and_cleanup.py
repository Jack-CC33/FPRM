"""Create metadata-bearing TIFF copies using private tag 65000.
Keep original files and sidecars untouched.
"""
from pathlib import Path
import json, re
import numpy as np
import tifffile

TAG = 65000

def distribute(folder):
    folder=Path(folder)
    def read(name):
        p=folder/name
        return json.loads(p.read_text(encoding='utf-8-sig')) if p.is_file() else None
    experiment=read('metadata.json')
    for image in sorted(folder.glob('I_*.tif')):
        b=image.stem[2:]
        base={'basename':b,'experiment':experiment,'cell_segmentation':read(f'FPRM_Step2_Metadata_{b}.json'),'mitochondrial_segmentation':read(f'FPRM_Step3_MitoPredictionMetadata_{b}.json'),'compartment_generation':read(f'FPRM_Step3_Metadata_{b}.json')}
        # Preserve original TIFFs; create metadata-bearing copies in a separate directory.
        out=folder/'MetadataEmbedded'; out.mkdir(exist_ok=True)
        def export(src, record):
            dst=out/src.name
            if dst.exists(): raise FileExistsError(dst)
            with tifffile.TiffFile(src) as t:
                if len(t.pages)!=1: raise ValueError('Expected single-page TIFF')
                p=t.pages[0]; arr=p.asarray()
                kwargs={'photometric':p.photometric,'metadata':None,'description':p.description or None}
                if p.colormap is not None: kwargs['colormap']=p.colormap
                if 282 in p.tags and 283 in p.tags:
                    kwargs['resolution']=(p.tags[282].value,p.tags[283].value)
                    if 296 in p.tags: kwargs['resolutionunit']=p.tags[296].value
                if p.extrasamples: kwargs['extrasamples']=p.extrasamples
                payload=json.dumps({'FPRM':record},ensure_ascii=True,allow_nan=False)
                tags=[(TAG,'s',0,payload,False)]
                for code in (50838,50839):
                    if code in p.tags:
                        tag=p.tags[code]; tags.append((code,tag.dtype,tag.count,tag.value,False))
                tifffile.imwrite(dst,arr,extratags=tags,**kwargs)
            if not np.array_equal(arr,tifffile.imread(dst),equal_nan=True):
                dst.unlink(); raise RuntimeError('Pixel verification failed')
        for prefix in ('FPRM_Step2_CellLabels_','FPRM_Step3_MitoPrediction_','FPRM_Step3_CellMask_','FPRM_Step3_MitoMask_','FPRM_Step3_CytoMask_','FPRM_Step3_SegmentationReview_'):
            p=folder/(prefix+b+'.tif')
            if p.is_file() and not (out/p.name).exists(): export(p,base)
        prefix='FPRM_Step4_AnalysisMetadata_'+b
        for p in sorted(folder.glob(prefix+'*.json')):
            suffix=p.stem[len(prefix):]
            if suffix and not re.fullmatch(r'_[2-9][0-9]*|_1[0-9]+',suffix): continue
            analysis=read(p.name)
            if analysis.get('basename')!=b or analysis.get('run_suffix')!=suffix: continue
            record={**base,'analysis':analysis}
            complete=folder/f'FPRM_CompleteMetadata_{b}{suffix}.json'
            if complete.exists(): print(f'Skip existing {complete.name}'); continue
            for kind,ext in [('CellID','.tif'),('SP','.tiff'),('FreeNADH','.tif'),('FreeNADH_Colormapped','.tif')]:
                src=folder/f'FPRM_Step4_{kind}_{b}{suffix}{ext}'
                if not src.is_file(): raise FileNotFoundError(src)
                export(src,record)
            complete.write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')

if __name__=='__main__':
    from tkinter import Tk,filedialog
    root=Tk(); root.withdraw()
    try: folder=filedialog.askdirectory(title='Select experiment folder')
    finally: root.destroy()
    if folder: distribute(folder)
