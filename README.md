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

Open the Step 4 macro in Fiji and run it with the measured calibration standards. Independent channel intensity support and the existing 16-bit conversion are retained pending further testing. CMRatio is relative free NADH cyto/mito; undefined ratios exclude the cell. Run `Tests/test_intensity_preservation.ijm` in Fiji with Bio-Formats and a representative two-channel CZI to compare conversion modes. The test changes the conversion preference to scaling disabled at completion. No CZI preparation automation is included yet.

This update is a merge overlay, not a complete repository: retain your existing license, general environment files, model and conversion utility. Remove the superseded notebook and mask-generation macro if present. No orchestration or CLI parser is included.
