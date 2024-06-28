

#@ File (label = "Input folder", style = "directory") input_dir
#@ File (label = "Output folder for segmentation", style = "directory") output_dir
#@ File (label = "LabKit classifier", style = "file") cilia_classifier


file_list = getFileList(input_dir)
number_of_files = file_list.length;
for (i=2; i<number_of_files; i++) {
	print("running file " + file_list[i]);
	close("*");
	open(input_dir + File.separator + file_list[i]);
	print("File loaded");
	title = getTitle();
	filename_only = substring(title , 0, lastIndexOf(title , '.'));
	//MIP
	run("Z Project...", "projection=[Max Intensity]");
	run("Split Channels");
	
	//Substract background
	selectWindow("C1-MAX_"+title);
	run("Subtract Background...", "rolling=200 sliding");
	
	//image min max normalization
	getMinAndMax(min, max);
	run("32-bit");
	run("Subtract...", "value="+min);
	run("Divide...", "value="+(max - min));

	//save(output_dir + File.separator + "normalized_" + filename_only + ".tif");
	print("pre-processing done, starting classification");
	//labkit stuff
	run("Segment Image With Labkit", "input=Composite segmenter_file=" + cilia_classifier + " use_gpu=true");
	rename("classified_image");
	
	print("Classification complete, saving files");
	run("Properties...", "channels=1 slices=1 frames=1 pixel_width=0.4151329 pixel_height=0.4151329 voxel_depth=0.9856603");
	save(output_dir + File.separator + "classified_" + filename_only + ".tif");
	close("*");
}

print("pipeline complete");