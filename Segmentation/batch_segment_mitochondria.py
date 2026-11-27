from pathlib import Path
from segment_mitochondria import segment
from generate_compartment_masks import generate_compartment_masks

def main():
    from tkinter import Tk, filedialog
    root = Tk(); root.withdraw()
    try:
        folder = filedialog.askdirectory(title='Select experiment folder')
        if not folder: return
        model = filedialog.askopenfilename(title='Select Ilastik model', filetypes=[('Ilastik','*.ilp')])
        if not model: return
    finally: root.destroy()
    folder = Path(folder); meta = folder/'metadata.json'; meta = meta if meta.is_file() else None
    for p in sorted(folder.glob('I_*.tif')):
        b = p.stem[2:]
        if not all((folder/n).is_file() for n in [f'FPRM_Step2_CellLabels_{b}.tif',f'FPRM_Step2_CellSegmentationOutlines_{b}_rois.zip',f'FPRM_Step2_Metadata_{b}.json']):
            print(f'{p.name}: missing Step 2 outputs'); continue
        try:
            segment(model,p,metadata_path=meta)
            generate_compartment_masks(p,metadata_path=meta)
        except Exception as e: print(f'{p.name}: skipped/failed: {e}')
if __name__ == '__main__': main()
