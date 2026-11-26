// Step 4: run after Python compartment generation. Independent channel support retained.
////////   User options    ////////////////////////////////
// for log file
date = getTime();
user = getInfo("user.name");

// options
// Limit between lower (vitrified) & upper (fast rotation) bounds.  Applied to Relative Free NADH values.
pixel_lim_opt = 1;											// 0 to map grayscale 0-255 between SPvit & SPfree, 1 to map to defined limits, specified below
User_pixel_lims = newArray(0.0, 0.6);						// If pixel_lim_option set to 1, uses this relative free NADH range for grayscale image colormap
function ColorMap(RelativeFreeNADH_Measured, option, LimLow, LimHigh){
	// If option for cmap limits is 0 (option==0), use [0.0, 1.0] as 0-255 gray value range
	// If option==1, use user-defined limits as 0-255 gray value range
	px_low = 0;
	px_high = 255;
	
	// Set display range
	Low = 0.0;
	High = 1.0;
	if(option) {
		Low = LimLow;
		High = LimHigh;
	}
	
	// Scale Relative Free NADH between limits
	fraction = RelativeFreeNADH_Measured;
	scaled_fraction = (fraction-Low)/(High-Low);
	
	// Confine measurement to bounds for colormapped image
	if(scaled_fraction < 0) {
		scaled_fraction = 0;
	}
	if(scaled_fraction > 1) {
		scaled_fraction = 1;
	}
	
	// Convert scaled Relative Free NADH value to 0-255 range
	PxVal = px_low + (px_high-px_low)*scaled_fraction;
	
	return round(PxVal)
}





////////   MAIN    ////////////////////////////////
run("Misc...", "divide=NaN run");     		// makes both x/0 and 0/0 = to NaN
run("Conversions...", " ");					// If 8/16/32-bit images are converted, no scaling of values.  For handling Masks.

Dialog.create("FPRM calibration standards");
Dialog.addNumber("SPvit", 0.725, 3);
Dialog.addNumber("SPfree", 0.975, 3);
Dialog.show();
SPvit = Dialog.getNumber(); SPfree = Dialog.getNumber();
if (isNaN(SPvit) || isNaN(SPfree) || SPvit<=0 || SPfree<=SPvit || SPfree>1e30) exit("Invalid standards");

// Organize File & Directory names
path = File.openDialog("Choose Intensity 'I' File:");  // Select I (intensity) file
dir = File.getParent(path);                            // pull out the folder name
File.setDefaultDir(path);

I_fname = File.getName(path);                          // file name
dotposition = lastIndexOf(I_fname,".");
basename = substring(I_fname, 2, dotposition);         // base name of image to be analyzed (prefix removed)

// Needed file names
P_fname = dir + File.separator + "P_" + basename + ".tif";
S_fname = dir + File.separator + "S_" + basename + ".tif";

cellrois_fname = dir + File.separator + "FPRM_Step2_CellSegmentationOutlines_" + basename + "_rois.zip";

cellmask_fname = dir + File.separator + "FPRM_Step3_CellMask_" + basename + ".tif";
cytomask_fname = dir + File.separator + "FPRM_Step3_CytoMask_" + basename + ".tif";
mitomask_fname = dir + File.separator + "FPRM_Step3_MitoMask_" + basename + ".tif";

// Reserve a common run suffix across all outputs.
runNumber=1;
do {
 suffix="";
 if (runNumber>1) suffix="_"+runNumber;
// Save file names
FPRM_CellID_img_fname = dir + File.separator + "FPRM_Step4_CellID_" + basename + suffix + ".tif";
FPRM_SP_img_fname = dir + File.separator + "FPRM_Step4_SP_" + basename + suffix + ".tiff";
FPRM_FreeNADH_img_fname = dir + File.separator + "FPRM_Step4_FreeNADH_" + basename + suffix + ".tif";
FPRM_FreeNADH_Colormapped_img_fname = dir + File.separator + "FPRM_Step4_FreeNADH_Colormapped_" + basename + suffix + ".tif";
FPRM_Results_fname = dir + File.separator + "FPRM_Step4_Results_" + basename + suffix + ".txt";
FPRM_Condensed_Results_fname = dir + File.separator + "FPRM_Step4_Condensed_Results_" + basename + suffix + ".txt";

 collision = File.exists(FPRM_CellID_img_fname) || File.exists(FPRM_SP_img_fname) || File.exists(FPRM_FreeNADH_img_fname) || File.exists(FPRM_FreeNADH_Colormapped_img_fname) || File.exists(FPRM_Results_fname) || File.exists(FPRM_Condensed_Results_fname);
 runNumber++;
} while (collision);
// Make record of how calculations were created (for recreation or debugging if necessary later)
info_str = "FPRM_Step4_SingleCellCalculations.ijm Macro\t" + date + "\t" + user + "\t" + basename;
File.append(info_str, FPRM_Results_fname);
File.saveString("", FPRM_Condensed_Results_fname);
info_str2 = "FPRM calibration standards: SPvit=" + SPvit + "; SPfree=" + SPfree;
File.append(info_str2, FPRM_Results_fname);
info_str3 = "           Color-coded image scaling from ";
if (pixel_lim_opt) {info_str3 = info_str3 + "provided limits:\t\t" + d2s(User_pixel_lims[0],2) + "\t" + d2s(User_pixel_lims[1],2);}
else {info_str3 = info_str3 + "standards' range:\t\t" + d2s(SPvit,3) + "\t" + d2s(SPfree,3);}
File.append(info_str3, FPRM_Results_fname);

///////////////////////////////////////////////////////////////////////////
// Open needed files & Create Masked Images
open(P_fname); // 32-bit
run("Select None");
rename("P");

open(S_fname); // 32-bit
run("Select None");
rename("S");

getDimensions(nx, ny, ch, slices, nz);

// cell mask
open(cellmask_fname); // 8-bit
run("Select None");
run("Divide...", "value=255.00000 ");
rename("Cell Mask");

roiManager("reset");
roiManager("Open", cellrois_fname);
n_roi = RoiManager.size;

run("Conversions...", "scale calibrate");
selectWindow("P");
run("16-bit");
selectWindow("S");
run("16-bit");
run("Conversions...", " ");

imageCalculator("Multiply create 16-bit", "P", "Cell Mask");
rename("P Cell Body");
imageCalculator("Multiply create 16-bit", "S", "Cell Mask");
rename("S Cell Body");


// cytoplasm mask
open(cytomask_fname); // 8-bit
run("Select None");
run("Divide...", "value=255.00000");
rename("Cyto Mask");

imageCalculator("Multiply create 16-bit", "P", "Cyto Mask");
rename("P Cyto");
imageCalculator("Multiply create 16-bit", "S", "Cyto Mask");
rename("S Cyto");

open(mitomask_fname); // 8-bit
run("Select None");
run("Divide...", "value=255.00000");
rename("Mito Mask");

imageCalculator("Multiply create 16-bit ", "P", "Mito Mask");
rename("P Mito");
imageCalculator("Multiply create 16-bit", "S", "Mito Mask");
rename("S Mito");





///////////////////////////////////////////////////////////////////////////
// Measure Full Cell average P, average S, P/S ratio, error of P/S, and P+S Integrated Intensity
run("Set Measurements...", "area mean centroid integrated limit redirect=None decimal=3");

newImage("SP_raw_values", "32-bit grayscale-mode", nx, ny, 1, 1, nz);
newImage("RelativeFreeNADH_raw_values", "32-bit grayscale-mode", nx, ny, 1, 1, nz);
newImage("RelativeFreeNADH_colormapped", "8-bit grayscale-mode", nx, ny, 1, 1, nz);

// iterating over multiple cells per FOV
for(roi = 0; roi < n_roi; roi++) {
	log_str = "ROI: " + d2s(roi,0);
	print(log_str);
	
	// FULL CELL MEASUREMENT
	selectWindow("P Cell Body");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Pcell = getResult("IntDen");
	Parea = getResult("Area");
	
	Xcentroid = getResult("X");
	Ycentroid = getResult("Y");
	run("Clear Results");
	
	selectWindow("S Cell Body");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Scell = getResult("IntDen");
	run("Clear Results");
	
	if (Parea<30 || Pcell<=0 || Scell<=0 || isNaN(Pcell) || isNaN(Scell)){
        File.append("Excluded ROI " + roi + ": invalid cell", FPRM_Results_fname);
		continue;
	}
	
	ValidP_cell = Parea;
    selectWindow("Cell Mask");
    roiManager("Select", roi);
    setThreshold(1,1,"raw");
    run("Measure");
    Area_cell = getResult("Area");
    resetThreshold();
    run("Clear Results");
	SP_cell = Scell / Pcell;
	RelativeFreeNADH_cell = (SP_cell - SPvit) / (SPfree - SPvit);
	Icell = Pcell + Scell;
	
	
	// CYTO MEASUREMENT
	selectWindow("P Cyto");
	run("Conversions...", "scale");
	run("16-bit");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Pcyto = getResult("IntDen");
	Parea = getResult("Area");
	run("Conversions...", " ");
	run("Clear Results");
	
	selectWindow("S Cyto");
	run("Conversions...", "scale");
	run("16-bit");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Scyto = getResult("IntDen");
	run("Conversions...", " ");
	run("Clear Results");
	
	if (Parea==0 || Pcyto<=0 || Scyto<=0 || isNaN(Pcyto) || isNaN(Scyto)){
        File.append("Excluded ROI " + roi + ": invalid cytoplasm", FPRM_Results_fname);
		continue;
	}
	
	ValidP_cyto = Parea;
    selectWindow("Cyto Mask");
    roiManager("Select", roi);
    setThreshold(1,1,"raw");
    run("Measure");
    Area_cyto = getResult("Area");
    resetThreshold();
    run("Clear Results");
	SP_cyto = Scyto / Pcyto;
	RelativeFreeNADH_cyto = (SP_cyto - SPvit) / (SPfree - SPvit);
	Icyto = Pcyto + Scyto;
	
	RelativeFreeNADH_cyto_color = ColorMap(RelativeFreeNADH_cyto, pixel_lim_opt, User_pixel_lims[0], User_pixel_lims[1]);
	
	// MITO MEASUREMENT
	selectWindow("P Mito");
	run("Conversions...", "scale");
	run("16-bit");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Pmito = getResult("IntDen");
	Parea = getResult("Area");
	run("Conversions...", " ");
	run("Clear Results");
	
	selectWindow("S Mito");
	run("Conversions...", "scale");
	run("16-bit");
	setThreshold(1, 65535, "raw");
	roiManager("Select", roi);
	run("Measure");
	Smito = getResult("IntDen");
	run("Conversions...", " ");
	run("Clear Results");
	
	if (Parea==0 || Pmito<=0 || Smito<=0 || isNaN(Pmito) || isNaN(Smito)){
		File.append("Excluded ROI " + roi + ": invalid mitochondria", FPRM_Results_fname);
		continue;
	}
	
	ValidP_mito = Parea;
    selectWindow("Mito Mask");
    roiManager("Select", roi);
    setThreshold(1,1,"raw");
    run("Measure");
    Area_mito = getResult("Area");
    resetThreshold();
    run("Clear Results");
	SP_mito = Smito / Pmito;
	RelativeFreeNADH_mito = (SP_mito - SPvit) / (SPfree - SPvit);
	Imito = Pmito + Smito;
	
	Pcell = Pcyto + Pmito;
	Scell = Scyto + Smito;
	
	RelativeFreeNADH_mito_color = ColorMap(RelativeFreeNADH_mito, pixel_lim_opt, User_pixel_lims[0], User_pixel_lims[1]);
	
    if (Area_cell<=0 || Area_mito<=0 || Area_cyto<=0 || RelativeFreeNADH_mito==0 || isNaN(RelativeFreeNADH_mito)) {
        File.append("Excluded ROI " + roi + ": undefined CMRatio", FPRM_Results_fname);
        continue;
    }
	// ROI - SPcyto - Cytoplasm intensity
	Cytoplasm_str = d2s(roi,0) + "\tCyto\tSPcyto\t" + d2s(SP_cyto,3) + "\tCyto Intensity\t" + d2s(Icyto,0) + "\tRelative free NADH\t" + d2s(RelativeFreeNADH_cyto,3) + "\tPixel Value\t" + d2s(RelativeFreeNADH_cyto_color,0);
	File.append(Cytoplasm_str, FPRM_Results_fname);
	
	
	
	// raw S/P value image
	selectWindow("Cell Mask");
	run("Duplicate...", "title=MaskCopy ignore");
	setBackgroundColor(0,0,0);
	roiManager("Select", roi);
	run("Clear Outside");
	
	imageCalculator("Multiply create 32-bit", "Cyto Mask", "MaskCopy");  // this gives wrong size image from mask generation macro sometimes
	rename("CytoROI");
	run("Multiply...", "value=SP_cyto");
	
	run("Concatenate...", "  title=SP_raw_values image1=SP_raw_values image2=CytoROI");
	selectWindow("SP_raw_values");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("SP_raw_values");
	close();
	selectWindow("SUM_SP_raw_values");
	rename("SP_raw_values");
	
	// raw RelativeFreeNADH value image
	imageCalculator("Multiply create 32-bit", "Cyto Mask", "MaskCopy");  // this gives wrong size image from mask generation macro sometimes
	rename("CytoROI");
	run("Multiply...", "value=RelativeFreeNADH_cyto");
	
	run("Concatenate...", "  title=RelativeFreeNADH_raw_values image1=RelativeFreeNADH_raw_values image2=CytoROI");
	selectWindow("RelativeFreeNADH_raw_values");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("RelativeFreeNADH_raw_values");
	close();
	selectWindow("SUM_RelativeFreeNADH_raw_values");
	rename("RelativeFreeNADH_raw_values");
	
	// color-mapped RelativeFreeNADH image
	imageCalculator("Multiply create", "Cyto Mask", "MaskCopy");
	rename("CytoROI");
	run("Multiply...", "value=RelativeFreeNADH_cyto_color");
	run("8-bit");
	
	run("Concatenate...", "  title=RelativeFreeNADH_colormapped image1=RelativeFreeNADH_colormapped image2=CytoROI");
	selectWindow("RelativeFreeNADH_colormapped");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	run("8-bit");
	selectWindow("RelativeFreeNADH_colormapped");
	close();
	selectWindow("SUM_RelativeFreeNADH_colormapped");
	rename("RelativeFreeNADH_colormapped");
	
	log_str = "\tCyto: " + d2s(RelativeFreeNADH_cyto,2) + "         " + d2s(RelativeFreeNADH_cyto_color,0);
	print(log_str);

	
	// ROI - PSmito - Mitochondria intensity
	Mitochondria_str = d2s(roi,0) + "\tMito\tSPmito\t" + d2s(SP_mito,3) + "\tMito Intensity\t" + d2s(Imito,0) + "\tRelative free NADH\t" + d2s(RelativeFreeNADH_mito,3) + "\tPixel Value\t" + d2s(RelativeFreeNADH_mito_color,0);
	File.append(Mitochondria_str, FPRM_Results_fname);
	
	// raw S/P value image
	imageCalculator("Multiply create 32-bit", "Mito Mask", "MaskCopy");
	rename("MitoROI");
	run("Multiply...", "value=SP_mito");
	
	wait(10);
	run("Concatenate...", "  title=SP_raw_values image1=SP_raw_values image2=MitoROI");
	selectWindow("SP_raw_values");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("SP_raw_values");
	close();
	selectWindow("SUM_SP_raw_values");
	rename("SP_raw_values");
	
	// raw RelativeFreeNADH image
	imageCalculator("Multiply create 32-bit", "Mito Mask", "MaskCopy");  // this gives wrong size image from mask generation macro sometimes
	rename("MitoROI");
	run("Multiply...", "value=RelativeFreeNADH_mito");
	
	wait(10);
	run("Concatenate...", "  title=RelativeFreeNADH_raw_values image1=RelativeFreeNADH_raw_values image2=MitoROI");
	selectWindow("RelativeFreeNADH_raw_values");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("RelativeFreeNADH_raw_values");
	close();
	selectWindow("SUM_RelativeFreeNADH_raw_values");
	rename("RelativeFreeNADH_raw_values");
	
	// color-mapped RelativeFreeNADH image
	imageCalculator("Multiply create", "Mito Mask", "MaskCopy");
	rename("MitoROI");
	run("Multiply...", "value=RelativeFreeNADH_mito_color");
	run("8-bit");
	
	wait(10);
	run("Concatenate...", "  title=RelativeFreeNADH_colormapped image1=RelativeFreeNADH_colormapped image2=MitoROI");
	selectWindow("RelativeFreeNADH_colormapped");
	wait(10);
	run("Z Project...", "projection=[Sum Slices]");
	run("8-bit");
	selectWindow("RelativeFreeNADH_colormapped");
	close();
	selectWindow("SUM_RelativeFreeNADH_colormapped");
	rename("RelativeFreeNADH_colormapped");
	
	selectWindow("MaskCopy");
	close(); 
	
	log_str = "\tMito: " + d2s(RelativeFreeNADH_mito,2) + "         " + d2s(RelativeFreeNADH_mito_color,0);
	print(log_str);
	
	
	
	// Cyto/Mito RATIO MEASUREMENT
	CMRatio = RelativeFreeNADH_cyto / RelativeFreeNADH_mito;

	// ROI - Total Integrated Intensity - Cellbody avg S/P - Relative Free NADH (cell avg) - SPcyto/SPmito
	Cellbody_str = d2s(roi,0) + "\tCell\tSPcell\t" + d2s(SP_cell,3) + "\tCell Intensity\t" + d2s(Icell,0) + "\tRelative free NADH\t" + d2s(RelativeFreeNADH_cell,3) + "\tCytofree/Mitofree\t" + d2s(CMRatio,3);
	File.append(Cellbody_str, FPRM_Results_fname);
    File.append("ROI " + roi + " mitochondrial area fraction=" + Area_mito/Area_cell + "; P support areas=" + ValidP_cell + "," + ValidP_mito + "," + ValidP_cyto, FPRM_Results_fname);
	
	// Condensed Results String
	// Cell # -- Cell Area -- Cell centroid X -- Y -- Cell Int -- Cell S/P -- Cell relative free NADH -- CytoMitoRatio (RelativeFreeNADH_cyto/RelativeFreeNADH_mito) -- Mito Area -- Mito Int -- Mito S/P -- Mito relative free NADH -- Cyto Area -- Cyto Int -- Cyto S/P -- Cyto relative free NADH
	Condensed_str = d2s(roi,0) + "\t" + d2s(Area_cell,0) + "\t" + d2s(Xcentroid,3) + "\t" + d2s(Ycentroid,3) + "\t" + d2s(Icell,0) + "\t" + d2s(SP_cell,3) + "\t" + d2s(RelativeFreeNADH_cell,3) + "\t" + d2s(CMRatio,3) + "\t";
	Condensed_str = Condensed_str + d2s(Area_mito,0) + "\t" + d2s(Imito,0) + "\t" + d2s(SP_mito,3) + "\t" + d2s(RelativeFreeNADH_mito,3) + "\t";
	Condensed_str = Condensed_str + d2s(Area_cyto,0) + "\t" + d2s(Icyto,0) + "\t" + d2s(SP_cyto,3) + "\t" + d2s(RelativeFreeNADH_cyto,3);
	File.append(Condensed_str, FPRM_Condensed_Results_fname);
}


// Save results
open(path);
run("Select None");
run("16-bit");
rename("I");
roiManager("Show All with labels");
run("Flatten");
selectWindow("I-1");
rename("Cell ID");
saveAs("Tiff", FPRM_CellID_img_fname);

selectWindow("SP_raw_values");
saveAs("Tiff", FPRM_SP_img_fname);

selectWindow("RelativeFreeNADH_raw_values");
saveAs("Tiff", FPRM_FreeNADH_img_fname);

selectWindow("RelativeFreeNADH_colormapped");
saveAs("Tiff", FPRM_FreeNADH_Colormapped_img_fname);


// Close all windows
waitForUser("Inspect / Save any images, then press OK to close all windows:");
run("Close All");