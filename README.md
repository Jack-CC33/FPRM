# FPRM
Analysis and utility resources for NAD(P)H fluorescence polarization ratiometric microscopy (FPRM)





## License and intended use

This software is provided to support scientific research, education, and the free exchange of knowledge within the academic community.

The software is made available under the **PolyForm Noncommercial License 1.0.0**. Researchers, educators, students, and other noncommercial users are encouraged to use, modify, adapt, and redistribute the code in accordance with the terms of this license.

Commercial use, including incorporation into commercial products or services, requires separate written permission from the copyright holder.

## Citation and acknowledgment

If you find this software useful in your research, teaching, or other academic activities, citation or acknowledgment of this repository is appreciated but is not a condition of use under the software license.

Users are encouraged to share improvements, extensions, and adaptations with the broader scientific community, although doing so is not required.

## AI-assisted development

Portions of this repository were developed with assistance from large language models (LLMs), including OpenAI's ChatGPT 6.2. These tools were used as collaborative aids in code development and documentation.

The human developers retain responsibility for the scientific concepts, implementation decisions, validation, and interpretation of results. The use of AI-assisted development tools does not imply endorsement by their respective developers or organizations.

## Disclaimer

This software is provided as is, without warranty. Users are responsible for independently verifying its suitability, accuracy, and reliability for their intended scientific or educational applications.

## Workflow

This method is intended for cells and tissue in which a user has already validated that ~ 440 nm autofluorescence signal highly correlates with a mitochondrial stain. See discussion of Fig 4 in https://doi.org/10.1126/sciadv.aed6844 for details. Analytical bias and mis-interpreted results are possible for non-validated samples.

Prepare aligned single-channel 2D I_, P_, S_ TIFFs with identical arbitrary basenames and .tif extensions. On the dedicated Cellpose computer activate your existing CP-SAM environment and run `python Segmentation/segment_cells.py`. Transfer the instance labels, ROI ZIP and JSON to the original image directory on the analysis computer. Instance labels are integers; binary masks are 0/255.

On Ubuntu create the environment with `conda env create -f CondaEnvironment/mito-segmentation.yaml`, activate `mito-segmentation`, and install Ilastik separately. Copy `ilastik_path.example.txt` to `ilastik_path.txt` and enter the actual launcher path. Run `python Segmentation/segment_mitochondria.py` and select the trained model and I image. Then run `python Segmentation/generate_compartment_masks.py` and select the same I image. No Cellpose installation is needed here. The review TIFF uses identical display scaling in both panels, yellow boundaries and cyan alpha 0.4. Metadata includes hashes; missing upstream metadata is recorded as null. Do not use training annotation TIFFs as predictions.

Open the Step 4 macro in Fiji and select Step 1 metadata.json to load calibration standards, or uncheck the option to enter them manually. Independent channel intensity support and the existing 16-bit conversion are retained pending further testing. CMRatio is relative free NADH cyto/mito; undefined ratios exclude the cell. Run `Tests/test_intensity_preservation.ijm` in Fiji with Bio-Formats and a representative two-channel CZI to compare conversion modes. The test changes the conversion preference to scaling disabled at completion. Run `Preparation/ImagePreparation.ijm` in Fiji with Bio-Formats to prepare CZI files and create experiment-level `metadata.json`.

This update is a merge overlay, not a complete repository: retain your existing license, general environment files, model and conversion utility. Remove the superseded notebook and mask-generation macro if present. No CLI parser is included.

### Optional experiment metadata

Each interactive Python segmentation script offers Step 1 `metadata.json` before image selection. Decline to use the original manual workflow. When selected, its current parent directory (not the recorded absolute path) initializes the I-image dialog. Select an I image in that directory. Copy metadata.json alongside the images when moving between Windows and Ubuntu. Programmatic callers may pass `metadata_path=...`; no automatic metadata discovery occurs. The segmentation records include the metadata filename and SHA-256 hash (null if omitted). Step 4 records the actual calibration values and their source. The ImageJ reader is limited to the numeric fields in the generated Step 1 schema, not arbitrary JSON.

Run `python -m unittest discover -s Tests -v` to test the Python metadata and mask logic.

### Batch segmentation and final metadata

Run `python Segmentation/batch_segment_cells.py` in the Cellpose environment, then `python Segmentation/batch_segment_mitochondria.py` in the mitochondrial environment. Select the experiment directory and model as prompted. Step 4 can be run interactively or with `Analysis/FPRM_BatchSingleCellCalculations.ijm`. Close existing image windows first. Batch loads calibration once, skips missing prerequisites, and stops on other macro errors. After analysis, run `python Segmentation/metadata_distribution_and_cleanup.py` in the mitochondrial environment to generate comprehensive JSON and metadata-bearing TIFF copies. Retain originals until TIFF calibration and display preservation have been validated.

### Single-cell results

`FPRM_Step4_Single_Cell_Results_basename.txt` is tab-delimited with the following 23 columns in order. `Cell_Index` is the existing zero-based ROI index; excluded cells can leave gaps.

1. `Cell_Index`
2. `Cell_Area`
3. `Cell_Centroid_X`
4. `Cell_Centroid_Y`
5. `Cell_Intensity`
6. `Cell_SP`
7. `Cell_RelativeFreeNADH`
8. `CytoMito_RelativeFreeNADH_Ratio`
9. `Mito_Area`
10. `Mito_Intensity`
11. `Mito_SP`
12. `Mito_RelativeFreeNADH`
13. `Cyto_Area`
14. `Cyto_Intensity`
15. `Cyto_SP`
16. `Cyto_RelativeFreeNADH`
17. `Mito_AreaFraction`
18. `Cell_ValidP_Area`
19. `Mito_ValidP_Area`
20. `Cyto_ValidP_Area`
21. `Mito_ColorValue`
22. `Cyto_ColorValue`

23. `Area_Units`

Areas use ImageJ spatial calibration: `px` for pixel counts, `um2` for square micrometers, or the spatial unit followed by `2`. P, S and masks must have matching spatial calibration; analysis stops on mismatch. The centroid is measured over P-positive cell pixels. Mito_AreaFraction is dimensionless. The minimum P area remains 30 in the current area units. Geometric areas retain fractional precision. Excluded cells and reasons are recorded in Step 4 AnalysisMetadata JSON. See [TIFF metadata guide](TIFF_METADATA.md).
