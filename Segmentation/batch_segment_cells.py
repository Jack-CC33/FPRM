from pathlib import Path
from segment_cells import segment_cells

def main():
    from tkinter import Tk, filedialog, messagebox
    root = Tk(); root.withdraw()
    try:
        folder = filedialog.askdirectory(title='Select experiment folder')
        if not folder: return
        model = None
        if not messagebox.askyesno('Cellpose', 'Use default CPSAM model?'):
            model = filedialog.askopenfilename(title='Select Cellpose model weights')
            if not model: return
        gpu = messagebox.askyesno('Cellpose', 'Use GPU?')
    finally: root.destroy()
    folder = Path(folder); meta = folder/'metadata.json'
    for p in sorted(folder.glob('I_*.tif')):
        try: segment_cells(p, metadata_path=meta if meta.is_file() else None, model_path=model, gpu=gpu)
        except Exception as e: print(f'{p.name}: skipped/failed: {e}')
if __name__ == '__main__': main()
