

print("pipeline starting")
from aicsimageio import AICSImage
import skimage
import os
from cellpose import models
import numpy as np
import tifffile
import pandas as pd
from tqdm import tqdm
import apoc
from scipy.spatial import distance
import yaml

print("imports complete")

CONFIG_NAME = 'config.yaml'

with open(CONFIG_NAME, "r") as f:
	config = yaml.safe_load(f)

scen_channel = config['scen_channel']
marker_channel = config['marker_channel']
nuclei_channel = config['nuclei_channel']
x_pad = config['x_pad']
y_pad = config['y_pad']
raw_data_folder = config['raw_data_folder']
classifications_folder = config['classifications_folder']
cl3D_filename = config['cl3D_filename']
output_folder = config['output_folder']

os.makedirs(output_folder, exist_ok=True)

#define functions

def ExtractSingleObject(image, centroid, x_pad = 50, y_pad = 50):
    cy = int(centroid[0])
    cx = int(centroid[1])

    subset_im = np.copy(image[cy-y_pad:cy+y_pad, cx-x_pad:cx+x_pad])
    return subset_im

def ExtractSingleObject3D(image, centroid, x_pad = 50, y_pad = 50):
    cy = int(centroid[0])
    cx = int(centroid[1])
    subset_im = np.copy(image[:, cy-y_pad:cy+y_pad, cx-x_pad:cx+x_pad])
    return subset_im

def filter_classified_image(im, lower_bound, upper_bound):
    '''takes a  label image, an upper and lower size bound and returns a labeled image conatining only the objects within that size bound'''
    props = skimage.measure.regionprops(im)
    for prop in tqdm(props, desc="Filtering senescent cell regions"):
        if prop['area'] < (lower_bound):
            im[im == prop['label']] = 0
        if prop['area'] > (upper_bound):
            im[im == prop['label']] = 0

    return im

def find_max_overlap(nuc_mask_subset, clean_scen_labs, min_overlap = 0.2):

    '''takes a celpose classified 3D nuclear mask image, and a scenescent label image with a single 
    label, and finds the corresponding nucleus based on maximum overlap
    returns nucleus label, nucleus centroid and overlap amount
    '''
    nuc_props = skimage.measure.regionprops(nuc_mask_subset)
    max_overlap = 0
    max_lab = - 1
    nuc_centroid = -1
    for nuc_prop in nuc_props:
        nuc_area = nuc_prop['area']
        overlap = np.count_nonzero(clean_scen_labs[nuc_mask_subset == nuc_prop['label']])
        if (overlap / nuc_area) > min_overlap:
            if overlap > max_overlap:
                max_overlap = overlap
                max_lab = nuc_prop['label']
                nuc_centroid = nuc_prop['centroid']
    return max_lab, nuc_centroid, max_overlap

def find_scenescent_cell_closest_to_center(scen_subset, min_area, x_pad=x_pad, y_pad=y_pad):
    '''takes a scenescent cell classified subset, a miniumum area, and the padding values and returns
    an image containing only the label closest to the center'''
    labs = skimage.measure.label(scen_subset)
    props = skimage.measure.regionprops(labs)
    min_dist = 9999
    min_prop = -1
    for prop in props:
        if prop['label'] > 1:       #do not include background lab
            if prop['area'] > min_area:
                centroid = (prop['centroid'][1], prop['centroid'][2])
                dist =  distance.euclidean(centroid, (x_pad,y_pad))
                if dist < min_dist:
                    min_dist = dist
                    min_prop = prop['label']
    if min_prop > 0:
        labs[labs != min_prop] = 0
    return labs

def classify_3D_crop(crop, erosion_num, cl3D_filename):
    '''takes a 3D crop and applies the 3D classifier
    dilates and erodes a set number of times (default = 4)
    returns the classified and dilated/eroded image'''

    clf = apoc.PixelClassifier(opencl_filename=cl3D_filename)
    result = clf.predict(image=crop)

    result_eroded = np.copy(result)
    for k in range(0, erosion_num):
        result_eroded = skimage.morphology.dilation(result_eroded)
    for k in range(0, erosion_num):
        result_eroded = skimage.morphology.erosion(result_eroded)

    return result_eroded

def find_matching_nucleus(nuc_centroid, nuc_crops_2D):
    '''Matches a 3D centroid of a nucleus to the closest nucleus in the 2D image'''
    nuc_props_2D = skimage.measure.regionprops(nuc_crops_2D)
    min_dist = 999999
    for prop in nuc_props_2D:
        dist =  distance.euclidean((nuc_centroid[1], nuc_centroid[2]), prop['centroid'])
        if dist < min_dist:
            min_dist = dist
            min_prop = prop['label']
    return min_prop, min_dist



  
model = models.Cellpose(model_type='cyto2')
files = os.listdir(raw_data_folder)

for filenum in range(0, len(files)):
    if filenum > -1:
        #load raw data and classified image
        print('loading images')
        filename = files[filenum][:-4]
        single_output_folder = output_folder + os.path.sep + filename
        os.makedirs(single_output_folder, exist_ok=True)
        print("Loading filename: " + filename)
   


        im = AICSImage(raw_data_folder + os.path.sep + filename + '.czi')
        raw_im = im.data[0,:,:,:,:]
        classified_im = skimage.io.imread(classifications_folder + os.path.sep + 'classified_' + filename + '.tif')

        #cellpose segment
        if os.path.isfile(single_output_folder + os.path.sep + filename + 'segmented_nuclei.tif'):
            print('segmentation image found, loading from disk')
            nuclei_im = skimage.io.imread(single_output_folder + os.path.sep + filename + 'segmented_nuclei.tif')
        else:
            print('starting 2D nuclear segmentation')
            im_to_segment = np.amax(raw_im[nuclei_channel, :, :], axis=0)
            nuclei_im, flows, styles, diams = model.eval(im_to_segment, diameter=20, flow_threshold=None)
            skimage.io.imsave(single_output_folder + os.path.sep + filename + '_segmented_nuclei.tif', nuclei_im, check_contrast=False)

        #erode and dilate the classified image
        eroded_lab_image = np.copy(classified_im)
        dilated = skimage.morphology.erosion(eroded_lab_image)
        for i in range(0, 4):
            eroded_lab_image = skimage.morphology.erosion(eroded_lab_image)
        for i in range(0, 4):
            eroded_lab_image = skimage.morphology.dilation(eroded_lab_image)

        #filter the classified objects

        nuc_labs = skimage.measure.label(nuclei_im)
        classified_labs = skimage.measure.label(eroded_lab_image)

        #calculate the minimum, mean, and maximum nuclear area
        areas = []
        all_nuc_props = skimage.measure.regionprops(nuc_labs)
        for nuc_prop in tqdm(all_nuc_props, desc="Measuring nuclei"):
            areas.append(nuc_prop['area'])
        mean_nuc_area = np.mean(areas)
        max_nuc_area = np.max(areas)
        min_nuc_area = np.min(areas)

        filtered_lab_image = filter_classified_image(classified_labs, mean_nuc_area * 2, max_nuc_area * 30)

        #extract crops
        crops = []
        nuc_crops = []
        centroids = []
        nuc_crops_2D = []

        scenescent_props = skimage.measure.regionprops(filtered_lab_image)
        for prop in tqdm(scenescent_props):
            crop_3D = ExtractSingleObject3D(raw_im[scen_channel,:,:,:], prop['centroid'], x_pad, y_pad)
            nuc_crop_3D = ExtractSingleObject3D(raw_im[nuclei_channel,:,:,:], prop['centroid'], x_pad, y_pad)
            nuc_crop_2D = ExtractSingleObject(nuclei_im, prop['centroid'], x_pad, y_pad)
            centroids.append(prop['centroid'])
            crops.append(crop_3D)
            nuc_crops.append(nuc_crop_3D)
            nuc_crops_2D.append(nuc_crop_2D)
        
        #loop through all cells in an image
        model = models.Cellpose(model_type='cyto2')
        matching_nuclei = np.zeros_like(nuclei_im)

        if os.path.isfile(single_output_folder  + os.path.sep + filename + "_senescent_nuclei.tif"):
            print('matched nuclei found, using saved image')
            matching_nuclei = skimage.io.imread(single_output_folder  + os.path.sep + filename + "_senescent_nuclei.tif")
        else:
            for cropnum, single_crop in enumerate(crops):
                os.system('clear')
                print("matching nuclei to cell " + str(cropnum+1) + " of " + str(len(crops)))
                #classify the 3D crop
                result_eroded = classify_3D_crop(single_crop, 4, cl3D_filename)
                #remove labs not near center
                labs = find_scenescent_cell_closest_to_center(result_eroded, 100, x_pad, y_pad)
                #cellpose segment the crops
                nuc_mask_subset, flows, styles, diams = model.eval(nuc_crops[cropnum], diameter=20, flow_threshold=None, do_3D=True)
                #find the nucleus with the maxium overlap
                max_lab, nuc_centroid, max_overlap = find_max_overlap(nuc_mask_subset, labs)
                if max_lab > 0:
                    #find the matching nucleus in the 2D classification
                    min_prop, min_dist = find_matching_nucleus(nuc_centroid, nuc_crops_2D[cropnum])
                    matching_nuclei[nuclei_im == min_prop] = min_prop
            print('saving matched nuclear image')
            skimage.io.imsave(single_output_folder  + os.path.sep + filename + "_senescent_nuclei.tif", matching_nuclei)
        #measure the nuclei
        matched_nuc_props = skimage.measure.regionprops(matching_nuclei)
        measurement_im = np.amax(raw_im[marker_channel,:,:,:], axis=0)

        intensities = []
        senescent = []
        labels = []

        for prop in tqdm(matched_nuc_props, desc="Quantifying intensity of senescent cells"):
            intensity = np.mean(measurement_im[matching_nuclei == prop['label']])
            intensities.append(round(intensity,2))
            labels.append(prop['label'])
            senescent.append(1)

        for prop in tqdm(all_nuc_props, desc="Quantifying intensity of non-senescent cells"):
            if prop['label'] not in labels:
                intensity = np.mean(measurement_im[nuc_labs == prop['label']])
                intensities.append(round(intensity,2))
                labels.append(prop['label'])
                senescent.append(0)

        output_df = pd.DataFrame()
        output_df['label'] = labels
        output_df['intensity'] = intensities
        output_df['senescent'] = senescent

        output_df.to_csv(single_output_folder + os.path.sep + filename + '_output.csv')


print("pipeline finished")