//  General note:
//  The Radial Profile Angle plugin is intended as a graphical interface plugin which displays all data automatically graphically.
//  Thus, it opens many windows during the measurements which is very memory consuming and slows down the measurements.
//  The additional windows can be supressed with setBatchMode(true), however this unfortunately creates a memory leak that I could not fix easily.



//  The image needs to be opened manually as a time stack, here of only one color channel.
//  Get image name and ID for later reference.
imageName = getTitle();
ide = getImageID();

//  Create a high contrast copy of the original stack to then convert into a binary and detect cells. 
//  This is only for better detection of the cells, not for the actual measurement which is done on the original image.
//  The values in setMinAndMax() may need to be adjusted depending on the properties of the image provided
run("Duplicate...", "duplicate");
setMinAndMax(0, 400);
setOption("ScaleConversions", true);
run("8-bit");
run("Smooth","stack");

//  Create the binary and use binary operations to first close holes in the black areas and then separate touching cells using the "Watershed" function.
//  For the images obtained here, one cycle of dilate worked well but the operations may need to be adjusted for other images.
setOption("BlackBackground", false);
run("Make Binary", "method=Otsu background=Default calculate create");
run("Dilate", "stack");
run("Fill Holes", "stack");
run("Watershed", "stack");

//  Here all measurements are activated which is not necessary but in case some of the information should be needed at some point afterwards, it is measured here and available.
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display add redirect=["+ imageName +"] decimal=9");

//  Detect the cells as particles based on their size and circularity. The values come from a smalll trial and error process of what is best suited to detect cells and may be adjusted. 
//  Detected particles that touch the border of the image are excluded and the overlay of all detected cells (ROIs) is displayed.
run("Analyze Particles...", "size=4.00-20.00 circularity=0.5-1.00 show=Overlay display exclude clear overlay add stack");

//  For each detected ROI the X and Y coordinates and the radius ("Major" is the diameter) are obtained and stored in new internal arrays to be assigned later to the radial intensity measurements 
n = nResults
x_coord =newArray(nResults)
y_coord = newArray(nResults)
radius = newArray(nResults)
for (i = 0; i < nResults(); i++) {
    x_coord[i] = getResult("X", i);
    y_coord[i] = getResult("Y", i);
    radius[i] = 0.5 * getResult("Major", i);
}

//  Save the results table to your favorite path and empty the results table afterwards.
saveAs("Results", "//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_1_eyfp/info.csv");
run("Clear Results");

//  Select the original image, i.e. the one that has not been changed in its contrast and smoothed to create the binary and apply the ROIs from the binary to it.
selectWindow(imageName);
run("From ROI Manager");

//  Get the pixel size for the image to convert the measured coordinates from µm to pixels, which is necessary for the "Radial Profile Angle" plugin.
//  w: width; h: height which have the same value obviously.
getPixelSize(micron, w, h);

//  Define two dummy indices to correctly arrange the results of the radial intensity measurements in the results table afterwards (s) 
//  and switch to the next slice of the stack once all cells on the current slice have been measured (b)
b = 1
s = 0

//  The following lines are the actual measurement of pixel intensity 

//  Loop through all detected ROIs present in the ROI manager 
for (k = 0; k < roiManager("count"); k++) {
	
	//  Select the k-th ROI
	roiManager("Select", k);
	
	//  Get current slice number
	a = getSliceNumber();
	
	//  If the current slice number a is larger than the previous slice number (b), the results are saved and the results table is cleared. 
	//  (s) is reset to 0 so the results of the following measurements are pasted starting with the first position of the results table again.
	//  This is done to cope with the high memory usage and the consequent slowing of the measurements that occurs over one slice of the stack.
	if (a > b) {
		saveAs("Results","//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_1_eyfp/radial-intensity/slice_"+b+".txt");
		run("Clear Results");
		s = 0;
	} 
	
	
	//  Use the "Radial Profile Angle" plugin, looping over 36 10° sectors to be measured. If the measured sector size is changed, the iterations can be adjusted here. 
	for (j = 0; j < 36; j++) {
		
		//  Execute the "Radial Profile Angle" plugin. the center as well as the radius of the cell to be measured are retrieved from the internal array created before (lines 38-46).
		//  The coordinates need to be converted from µm to pixels using the before defined width and height.
		//  The starting angle of the measured sector is shifted with the dummy j.
		//  Integration angle is set to 5° which is interpreted by the plugin as 5° clockwise and 5° counter-clockwise around the central starting angle leading in total to a 10° sector.
		//  The radius is rounded because the plugin often creates NaN at the maximum radii if they are not integer. 
		//  NaNs are also produced for very small radii but for this aplication they are filtered in the subsequent R analysis which focuses on the large radii that contain the membrane area of the cells.
		run("Radial Profile Angle", " x_center="+(1/w)*x_coord[k]+" y_center="+(1/h)*y_coord[k]+" radius="+round((1/h)*radius[k])+" starting_angle="+j * 10+" integration_angle=5 ");
		
		//  Use the plugin provided macro extension functions to extract values from the plugin measurements and paste them to the results table. 
		//  A summary of the Ext. functions can be found in the macro documentation at http://questpharma.u-strasbg.fr/html/radial-profile-ext.html
		//  setResult("column", line, value) is used to paste the correct values in the correct positions in the table.
		//  The exact number of loops here depends on the size of the curve array created by the plugin (Ext.getBinSize), i.e. the size of the cell that is measured.
		for (i = 0; i != Ext.getBinSize; i++) {
			setResult("angle", s + j * parseInt(Ext.getBinSize) + i, j * 10);
			setResult("cell", s + j * parseInt(Ext.getBinSize) + i, k + 1);
			setResult("sector.start", s + j * parseInt(Ext.getBinSize) + i, (j * 10) - 5);
			setResult("sector.end", s + j * parseInt(Ext.getBinSize) + i, (j * 10) + 5);
			setResult("radius", s + j * parseInt(Ext.getBinSize) + i, parseFloat(Ext.getXValue(i)));
			setResult("int.intensity.norm", s + j * parseInt(Ext.getBinSize) + i, parseFloat(Ext.getYValue(0, i)));
			setResult("timepoint", s + j * parseInt(Ext.getBinSize) + i, a - 1);
			setResult("X", s + j * parseInt(Ext.getBinSize) + i, x_coord[k]);
			setResult("Y", s + j * parseInt(Ext.getBinSize) + i, y_coord[k]); 
		}
		
		updateResults();
		
		//  The original image is selected as the active image and all other windows are closed. 
		//  This is because the plugin opens a graphical display of the measured data after each measurement which needs to be closed to conserve memory space.
		selectWindow(imageName);
		close("\\Others");
		
		//  Deselect the ROI which has been measured. 
		run("Select None");
	}
	//  Update the sorting dummies according to the current number of measurements and slice number in the stack after one cell has been measured completely.
	s += parseInt(Ext.getBinSize)*j;
	b = getSliceNumber();
	
	//  Save the results table in the case that the last ROI on the last slice has been measured.
	if (k == roiManager("count")-1) {
		saveAs("Results","//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_1_eyfp/radial-intensity/slice_"+b+".txt");
	
	//  Return to the  beginning of the loop and select the next ROI from the ROI manager.
	}
	
}
