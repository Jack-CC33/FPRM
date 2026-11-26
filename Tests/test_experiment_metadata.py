import hashlib
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import tifffile
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'Segmentation'))

class TestMetadata(unittest.TestCase):
    def test_loaders(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'metadata.json'
            p.write_text(json.dumps({'SPvit': .725, 'SPfree': .975, 'processing_directory': 'D:/old/location'}))
            for name in ('segment_cells', 'segment_mitochondria', 'generate_compartment_masks'):
                load = importlib.import_module(name).load_experiment_metadata
                self.assertEqual(load(None), (None, None))
                folder, ref = load(p)
                self.assertEqual(folder, p.parent.resolve())
                self.assertEqual(ref['sha256'], hashlib.sha256(p.read_bytes()).hexdigest())
            for values in ((0,1), (1,1), (1,.5), (float('nan'),1), ('0.725',.975)):
                p.write_text(json.dumps(dict(zip(('SPvit','SPfree'),values))))
                for name in ('segment_cells','segment_mitochondria','generate_compartment_masks'):
                    with self.assertRaises(ValueError):
                        importlib.import_module(name).load_experiment_metadata(p)

    def test_compartment_reference(self):
        from generate_compartment_masks import generate_compartment_masks
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)
            meta = folder / 'metadata.json'
            meta.write_text('{"SPvit":0.725,"SPfree":0.975}')
            p = folder / 'I_example.tif'
            l = folder / 'labels.tif'
            m = folder / 'prediction.tif'
            for f, a in ((p, [[10,20],[30,40]]),(l,[[0,1],[1,1]]),(m,[[0,255],[0,0]])):
                tifffile.imwrite(f,np.array(a,dtype=np.uint16))
            generate_compartment_masks(p,l,m,metadata_path=meta)
            result = json.loads((folder/'FPRM_Step3_Metadata_example.json').read_text())
            self.assertEqual(result['experiment_metadata']['sha256'],hashlib.sha256(meta.read_bytes()).hexdigest())
            other = folder/'other'; other.mkdir()
            with self.assertRaises(ValueError):
                generate_compartment_masks(other/'I_example.tif',l,m,metadata_path=meta)
