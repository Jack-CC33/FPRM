import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import tifffile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'Segmentation'))
from generate_compartment_masks import generate_compartment_masks

class TestMasks(unittest.TestCase):
    def test_partition(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'I_example.tif'; l=Path(d)/'labels.tif'; m=Path(d)/'prediction.tif'
            labels=np.array([[0,1,1],[0,2,2]],dtype=np.uint16)
            pred=np.array([[255,255,0],[255,0,0]],dtype=np.uint8)
            for f,x in ((p,labels*100),(l,labels),(m,pred)): tifffile.imwrite(f,x)
            out=generate_compartment_masks(p,l,m)
            c,mi,cy=[tifffile.imread(out[k]) for k in ('CellMask','MitoMask','CytoMask')]
            np.testing.assert_array_equal(c,(labels>0)*255)
            np.testing.assert_array_equal(mi,[[0,255,0],[0,0,0]])
            np.testing.assert_array_equal(c,mi|cy)
            self.assertFalse(np.any(mi&cy))
            self.assertEqual(tifffile.imread(out['SegmentationReview']).shape,(2,6,3))
            with self.assertRaises(FileExistsError): generate_compartment_masks(p,l,m)
            tifffile.imwrite(m,np.full(labels.shape,2,dtype=np.uint8))
            with self.assertRaises(ValueError): generate_compartment_masks(p,l,m,overwrite=True)
if __name__=='__main__': unittest.main()
