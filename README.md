<p align="center">
  <a href="https://ai4life.eurobioimaging.eu/open-calls/">
    <img src="https://github.com/ai4life-opencalls/.github/blob/main/AI4Life_banner_giraffe_nodes_OC.png?raw=true">
  </a>
</p>


# Project 47: Automated Detection and Quantification of Senescent Cells

A collaboration with Ana Filipa Pombinho Isidroe at Universidade de Lisboa. 

## Introduction

In this project, cross sections of mouse neural tube are stained with an antibody against a marker of senescence, as well as a nuclear marker that needs to be quantified in senescent and non-senescent cells. These input data are stitched 3D stacks of the entire cross-section. The primary difficulty is the inconsistency of the staining of the senescence marker, and the fact that there a protrusions emanating from the senescent cells which could overlap with other (non-senescent) nuclei, potentially confusing the analysis. 

In order to perform this analysis, we first flattened the 3D stack by sum projecting along the Z-axis. We then segmented all nuclei using [cellpose](https://cellpose.readthedocs.io/en/latest/api.html) and used the senescence channel to perform a first-pass pixel classification using [labkit](https://imagej.net/plugins/labkit/) and a classfier trained to detect senescent cell bodies, while excluding (as much as possible) protrusions.

Following this, for each detected senescent cell body (larger than a certain threshold), we extracted a 3D volume around the cell from the raw data and did a second (higher quality) 3D pixel classification in order to more precisely segment the senescent cell. From this, we are able to detect the nucleus with the maximal overlay with this mask in 3D. This 3D nucleus was then mapped back to the 2D nuclear segmented image by finding the closest nucleus between the 3D and 2D images.

From this 2D nucleus, the mean intensity of the marker channel was quantified, along with the non-senescent nuclei. These files are saved into a .csv file.

## Running the pipeline

### Part 1: Labkit classification

Run the script `labkit_classification.ijm` in Fiji. This script requires Fiji and Labkit to be installed. You can download Fiji from [here](https://imagej.net/software/fiji/downloads). If you have an old Fiji installation, [Labkit](https://imagej.net/plugins/labkit/) might not be there by default. To install it, go to `Help > Update` in Fiji. Then, click on `Manage Update Sites`. From that list, scroll down to `Labkit`, select the checkbox, and click the `Apply and Close` button. Then from the Updated, select `Apply changes`.

This script takes three input parameters:

`input folder:` The input folder where the raw data is stored

`output folder for segmentation:` The output data to store the segmentation ma

`labkit classifier:` The location of the classifier model for Labkit segmentation (by default, you should use `senescence_classifier.classifier` in the models folder)

The output of this script is, for each image, a semantic segmentation containing the masks for all detected senescent cells (as well as, possibly, some smaller detections that will be filtered out in subsequent steps)

### Step 2: Detection, matching and quantification of nuclei

Make a conda environment using the included yaml source file by typing `conda create -f environment.yaml`

Activate this environment by typing `conda activate senescence` 

If necessary, update the `config.yaml` file. In this file, you can set the following parameters

`scen_channel:` The position in the stack of the senescence marker

`marker_channel:` The position in the stack of the nuclear marker to be measured

`nuceli_channel:` The position in the stack of the DAPI staining

`x_pad:` How many pixels to extract in each direction of the x dimension around each detected 2D cell for 3D classification

`y_pad:` How many pixels to extract in each direction of the y dimension around each detected 2D cell for 3D classification

`raw_data_folder:` The location of the raw data (same folder as in step 1)

`classifications_folder:` The location of the output from step 1

`cl3D_filename:` The classifier file for 3D classification (default is `scen_3D.cl` in the `models` folder)

`output_folder:` The folder to store the output of the script

Next, run the script by typing `python3 analyze_nuclei.py`.

For each file in the input folder, the script will output a folder containing the following:

`<imagename>_segmented_nuclei.tif:` A tiff file containing the masks of all the segmented nuclei in 2D

`<imagename_senescent_nuclei.tif:` A tiff file containing only the masks of the nuclei which were matched to senescent cells, for validation purposes

`<imagename>_output.csv:` A csv file containing the mean intensities of all nuclei in the `marker_channel`. The csv has the following fields:

`label:` The label of the nuclei (corresponding to the value of the mask in the output tiff files

`intensity:` The mean intensity in the marker channel

`senescent:` A value of `1` indicates a nuclei that was matched to a senescent cell, while a value of `0` indicates a nucleus that was not matched to a senescent cell
