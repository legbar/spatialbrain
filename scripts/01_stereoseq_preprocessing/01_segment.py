import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
from scipy import ndimage as ndi

from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage import feature, measure, segmentation, color, exposure, filters, morphology, restoration, img_as_ubyte
from skimage.util import map_array

from PIL import Image as im

import pandas as pd

import glob

import spateo as st

from joblib import Parallel, delayed
from functools import partial

from cupyx.scipy import ndimage as cu_ndi
import cupy as cp

import os
os.chdir('/f_active/paper_23/') # laune
#os.chdir('/home/ubuntu/projects/paper_23/') # aws

def process_layer(adata, layer):
    if layer == 'X':
        m = adata.X
    else:
        m = adata.layers[layer]
    m = m.A
    print(np.sum(m))
    m = exposure.rescale_intensity(m)
    print('gaussian...')
    m = cp.array(m)
    m = cu_ndi.gaussian_filter(m, sigma=5)
    print('top hat...')
    m = cu_ndi.white_tophat(m, footprint = morphology.disk(50))
    m = cp.asnumpy(m)
    m = img_as_ubyte(exposure.rescale_intensity(m))
    print(m.dtype)
    return(m)

gems = {x.split('/')[-1].split('.tissue')[0]:x for x in glob.glob('input/01_stereoseq_preprocessing/gem/**/**/*gem.gz', recursive=True)}
spateo_in_files = {x.split('/')[-1].split('.counts')[0]:x for x in glob.glob('input/01_stereoseq_preprocessing/gem/**/**/*spateo.csv.gz', recursive=True)}

def segment(sn):

    try:
        adata = st.io.read_bgi_agg(spateo_in_files[sn])
        # [10000:11000,10000:11000]
    
        # gaussian filtering of each layer
        unspliced = process_layer(adata, 'unspliced')
        # spliced = process_layer(sn, 'spliced')
        X = process_layer(adata, 'X')
    
        print('otsu unspliced')
        thresh = filters.threshold_otsu(unspliced)
        binary_unspliced = unspliced >= thresh
    
        # thresh = filters.threshold_otsu(spliced)
        # binary_spliced = spliced >= thresh
    
        print('otsu X')
        thresh = filters.threshold_otsu(X)
        binary_X = X >= thresh
    
        # watershed on unspliced layer
        print('distance...')
        distance = ndi.distance_transform_edt(binary_unspliced)
        local_max_coords = feature.peak_local_max(distance, 
                                                  min_distance=7, 
                                                  labels = binary_unspliced)
        local_max_mask = np.zeros(distance.shape, dtype=bool)
        local_max_mask[tuple(local_max_coords.T)] = True
        print('label...')
        markers = measure.label(local_max_mask)
        print('watershed...')
        segmented_cells = segmentation.watershed(-distance, markers, watershed_line=True)
        print('dilation...')
        watershed_lines = morphology.binary_dilation(segmented_cells == 0, footprint=np.ones((5,5)))
    
        # generate post-watershed mask and fill holes
        mask = (binary_X + binary_unspliced).astype('int')
        mask[watershed_lines] = 0
        seed = np.copy(mask)
        seed[1:-1, 1:-1] = mask.max()
        print('erosion...')
        mask = morphology.reconstruction(seed, mask, method = 'erosion') # fill holes
        
        # filter on area
        print('measure...')
        label = measure.label(mask)
        label_overlay = color.label2rgb(label, image=mask, bg_label=0)
        df = pd.DataFrame(measure.regionprops_table(label, X, properties=['label', 'area', 'intensity_mean', 'equivalent_diameter_area', 'eccentricity', 'centroid']))
        keep = df['label'] * (df['equivalent_diameter_area'] > 15)
        print('map array...')
        label_filtered = map_array(label, np.asarray(df['label']), np.asarray(keep))
        label_filtered[label_filtered > 0] = 1
        mask = label_filtered
        
        # watershed on X
        print('distance...')
        distance = ndi.distance_transform_edt(mask)
        local_max_coords = feature.peak_local_max(distance, 
                                                  min_distance=7, 
                                                  labels = mask)
        local_max_mask = np.zeros(distance.shape, dtype=bool)
        local_max_mask[tuple(local_max_coords.T)] = True
        markers = measure.label(local_max_mask)
        print('watershed...')
        segmented_cells = segmentation.watershed(-distance, markers, watershed_line=True)
        print('dilation...')
        watershed_lines = morphology.binary_dilation(segmented_cells == 0, footprint=np.ones((3,3)))
        
        # remove X-watershed lines from mask
        mask[watershed_lines] = 0
        
        # filter on area
        print('measure...')
        label = measure.label(mask)
        label_overlay = color.label2rgb(label, image=mask, bg_label=0)
        df = pd.DataFrame(measure.regionprops_table(label, properties=['label', 'equivalent_diameter_area']))
        keep = df['label'] * (df['equivalent_diameter_area'] > 15)
        print('map array...')
        label_filtered = map_array(label, np.asarray(df['label']), np.asarray(keep))
        label_filtered[label_filtered > 0] = 1
        mask = label_filtered
        
        # filter on overlap with unspliced mask and eccentricity
        df = pd.DataFrame(measure.regionprops_table(label, binary_unspliced, properties=['label', 'intensity_max', 'eccentricity']))
        keep = df['label'] * ((df['intensity_max'] > 0) & (df['eccentricity'] < 0.8))
        label_filtered = map_array(label, np.asarray(df['label']), np.asarray(keep))
        label_filtered[label_filtered > 0] = 1
        mask = label_filtered
        
        print('save data...')
        data = im.fromarray((mask * 255).astype('uint8'))
        data.save('input/01_stereoseq_preprocessing/segmentation/masks/' + sn + '_mask.png')

        # edited 230928
        # data = im.fromarray(img_as_ubyte(exposure.rescale_intensity(unspliced)))
        # data.save('input/01_stereoseq_preprocessing/segmentation/masks/' + sn + '_unspliced.png')
    
        # data = im.fromarray(img_as_ubyte(exposure.rescale_intensity(X)))
        # data.save('input/01_stereoseq_preprocessing/segmentation/masks/' + sn + '_X.png')
    
        np.save('input/01_stereoseq_preprocessing/segmentation/masks/' + sn + '_mask.npy', mask)
    except:
      print(sn + " failed")

for sn in gems.keys():
    print(sn)
    segment(sn)

#Parallel(n_jobs=2)(delayed(segment)(sn) for sn in to_do)
