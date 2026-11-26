// Fiji + Bio-Formats required. Run in a clean Fiji session.
// Each Bio-Formats series is treated as one FOV (roi index).
files = newArray(0);
more = true;
while (more) {
    p = File.openDialog("Select CZI file");
    if (!endsWith(toLowerCase(p), ".czi")) { exit("Expected a CZI file"); }
    for (i=0; i<files.length; i++) { if (files[i]==p) { exit("Duplicate input file"); } }
    files = Array.concat(files, newArray(p));
    Dialog.create("Input files");
    Dialog.addCheckbox("Select another CZI file", false);
    Dialog.show(); more = Dialog.getCheckbox();
}
out = getDirectory("Choose or create processing directory");
if (File.exists(out+"metadata.json")) { exit("metadata.json already exists. Choose a new directory."); }
run("Bio-Formats Macro Extensions");
Ext.setId(files[0]); Ext.getSeriesCount(nseries); Ext.setSeries(0); Ext.getSizeC(nc); Ext.close();
if (nc<2) { exit("At least two channels required"); }
// Open first series for visual channel identification.
openSeries(files[0], 0);
preview = getImageID();
waitForUser("Inspect channels", "Use the channel slider to identify P and S, then click OK.");
choices = newArray(nc);
for (i=0; i<nc; i++) { choices[i] = ""+(i+1); }
Dialog.create("Channel mapping");
Dialog.addMessage("Inspect the first image using the channel slider before confirming.\nChannel numbers are 1-based.");
Dialog.addChoice("P channel", choices, "1");
Dialog.addChoice("S channel", choices, "2");
Dialog.show(); pc=parseInt(Dialog.getChoice()); sc=parseInt(Dialog.getChoice());
if (pc==sc) { exit("P and S must be distinct"); }
names = newArray(nc);
if (nc>2) {
    Dialog.create("Additional channel names");
    for (i=0; i<nc; i++) { if (i+1!=pc && i+1!=sc) { Dialog.addString("Channel "+(i+1), "ch"+(i+1)); } }
    Dialog.show();
    for (i=0; i<nc; i++) { if (i+1!=pc && i+1!=sc) { names[i]=Dialog.getString(); } }
}
names[pc-1]="P"; names[sc-1]="S";
for (i=0; i<nc; i++) {
    if (!matches(names[i], "[A-Za-z0-9][A-Za-z0-9_-]*")) { exit("Invalid channel name: "+names[i]); }
    if (toLowerCase(names[i])=="i") { exit("I is reserved for the sum"); }
    for (j=0; j<i; j++) { if (toLowerCase(names[i])==toLowerCase(names[j])) { exit("Duplicate channel name"); } }
}
selectImage(preview); close();
Dialog.create("Experiment metadata");
Dialog.addNumber("SPvit", 0.725, 3, 1, "");
Dialog.addNumber("SPfree", 0.975, 3, 1, "");
Dialog.addString("Biological condition / sample note", "", 60);
Dialog.addCheckbox("Use 32-bit sum (avoid integer saturation)", false);
Dialog.show(); vit=Dialog.getNumber(); free=Dialog.getNumber(); note=Dialog.getString(); floatSum=Dialog.getCheckbox();
if (isNaN(vit)||isNaN(free)||vit<=0||free<=vit||free>1e30) { exit("Invalid standards"); }
// Preflight every output name and channel count before writing anything.
planned=newArray(0);
for (f=0; f<files.length; f++) {
    Ext.setId(files[f]); Ext.getSeriesCount(ns);
    for (r=0; r<ns; r++) {
        Ext.setSeries(r); Ext.getSizeC(cc); Ext.getSizeZ(nz); Ext.getSizeT(nt);
        if (cc!=nc) { exit("Channel count differs: "+files[f]+" series "+r); }
        for (t=0; t<nt; t++) {
        for (z=0; z<nz; z++) {
            b=planeBase(files[f],t,z,r,nt,nz,ns);
            for (c=0; c<=nc; c++) {
                prefix="I"; if (c<nc) { prefix=names[c]; }
                target=out+prefix+"_"+b+".tif";
                if (File.exists(target)) { exit("Output exists: "+target); }
                for (k=0; k<planned.length; k++) { if (toLowerCase(planned[k])==toLowerCase(target)) { exit("Output name collision: "+target); } }
                planned=Array.concat(planned,newArray(target));
            }
        }
        }
    }
    Ext.close();
}
setBatchMode(true);
for (f=0; f<files.length; f++) {
    Ext.setId(files[f]); Ext.getSeriesCount(ns); Ext.close();
    for (r=0; r<ns; r++) {
        openSeries(files[f],r); src=getImageID(); getDimensions(w,h,cc,nz,nt);
        if (cc!=nc) { exit("Imported channel count differs from mapping"); }
        for (t=0; t<nt; t++) {
        for (z=0; z<nz; z++) {
            b=planeBase(files[f],t,z,r,nt,nz,ns);
            for (c=1; c<=nc; c++) {
                selectImage(src); Stack.setPosition(c,z+1,t+1);
                run("Duplicate...", "title=FPRM_plane duplicate channels="+c+" slices="+(z+1)+" frames="+(t+1));
                getDimensions(dw,dh,dc,dz,dt);
                if (dw!=w||dh!=h||dc!=1||dz!=1||dt!=1) { exit("Plane extraction failed"); }
                saveAs("Tiff",out+names[c-1]+"_"+b+".tif");
                if (c==pc) { pid=getImageID(); } else { if (c==sc) { sid=getImageID(); } else { close(); } }
            }
            operation="Add create"; if (floatSum) { operation="Add create 32-bit"; }
            imageCalculator(operation,pid,sid);
            saveAs("Tiff",out+"I_"+b+".tif"); close();
            selectImage(pid); close(); selectImage(sid); close();
        }
        }
        selectImage(src); close();
    }
}
setBatchMode(false);

json="{\n  \"schema_version\": \"1.0\",\n  \"processing_directory\": "+quote(out)+",\n  \"SPvit\": "+vit+",\n  \"SPfree\": "+free+",\n  \"channel_names\": {";
for (i=0; i<nc; i++) { 
	if (i>0) {
		json=json+", ";
		}
	json=json+quote(""+(i+1))+": "+quote(names[i]);
	}
json=json+"},\n  \"sample_note\": "+quote(note)+",\n  \"sum_mode\": ";
if (floatSum) {
	json=json+quote("32-bit"); 
	}
else {
	json=json+quote("native");
	}

json=json+",\n  \"source_files\": [";
for (i=0; i<files.length; i++) { 
	if (i>0) {
		json=json+", ";
		}
	json=json+quote(files[i]);
	}
json=json+"]\n}\n";
File.saveString(json,out+"metadata.json");
print("FPRM preparation complete: "+planned.length+" TIFFs in "+out);

function openSeries(path,r) {
    run("Bio-Formats Importer", "open=["+path+"] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+(r+1));
}
function planeBase(path,t,z,r,nt,nz,ns) {
    name=File.getName(path); b=substring(name,0,lastIndexOf(name,"."));
    if (nt>1) { b=b+"_t"+t; }
    if (nz>1) { b=b+"_z"+z; }
    if (ns>1) { b=b+"_roi"+r; }
    return b;
}
function quote(s) {
    result="\"";
    for (q=0; q<lengthOf(s); q++) {
        ch=substring(s,q,q+1); code=charCodeAt(s,q);
        if (ch=="\\" || ch=="\"") { result=result+"\\"+ch; } else {
            if (code<32) { hex=toHex(code); if (lengthOf(hex)==1) { hex="0"+hex; } result=result+"\\u00"+hex; }
            else { result=result+ch; }
        }
    }
    return result+"\"";
}
