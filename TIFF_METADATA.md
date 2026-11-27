# Reading FPRM TIFF metadata

Run `python Segmentation/metadata_distribution_and_cleanup.py` after analysis.
Select the experiment folder. Original TIFFs and sidecars remain untouched;
metadata-bearing copies are in `MetadataEmbedded`. Private TIFF tag 65000
contains ASCII-escaped JSON on the first page. It is not a standard ImageJ
metadata field and may be discarded by other software when saving images.

```python
import json
import tifffile

with tifffile.TiffFile('MetadataEmbedded/FPRM_Step3_MitoMask_example.tif') as tif:
    tag = tif.pages[0].tags.get(65000)
    if tag is None:
        raise ValueError('No embedded FPRM metadata')
    record = json.loads(tag.value)['FPRM']
    pixels = tif.asarray()

print(json.dumps(record, indent=2))
mito = record.get('mitochondrial_segmentation') or {}
print(mito.get('model'), mito.get('model_path'))
```

The record contains experiment, cell_segmentation, mitochondrial_segmentation,
compartment_generation and (only for Step 4 outputs) analysis. Missing sidecars
are represented by null. Analysis includes actual calibration, exclusions and
run suffix. The matching `FPRM_CompleteMetadata_<basename><suffix>.json` contains
all available records. Original file hashes refer to original files, not copies.
Do not regenerate segmentation between analysis and distribution: the script
cannot reconstruct historical segmentation provenance from overwritten sidecars.
Keep original files; verify calibration, LUTs and display in Fiji before relying
on copies. Existing embedded files are not overwritten.
