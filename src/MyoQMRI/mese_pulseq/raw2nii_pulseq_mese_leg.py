"""
    This code reconstructs the raw data which was acquired with the pulseq 
    MESE sequence.
    Data is reconstructed as [x,y,z,echoes] and saved as nifti.
    The seq file came together with a json file. Both must be in the same folder
    as the acquired data to run the reconstruction.
    Optional input parameters:
    -p --path       path to the dataset (default: current working directory)
"""

import os
import re
import pathlib
import numpy as np
from pypulseq import Sequence
import nibabel as nib
import json
from argparse import ArgumentParser
import mapvbvd
import twixtools as tx

### FUNCTIONS ###
# Read Siemens raw data in twix format
def read_raw_data(filename):
    # Read twix file
    twix_obj = mapvbvd.mapVBVD(filename)
    
    # If twix obj contains multiple scans, pick last one
    if isinstance(twix_obj, list):
        twix_obj = twix_obj[-1]
    
    # Load unsorted data: Do not remove oversampling and do not reflect lines according to the REV label
    twix_obj.image.removeOS = False
    twix_obj.image.disableReflect = True
    kdata = twix_obj.image.unsorted() # Shape: [N_adc, N_coils, N_meas]
    
    # Load in phasecor data and sort it into kdata by memory position (i.e. acquisition number)
    if hasattr(twix_obj, 'phasecor'):
        twix_obj.phasecor.removeOS = False
        twix_obj.phasecor.disableReflect = True
        kdata_phasecor = twix_obj.phasecor.unsorted() # Shape: [N_adc, N_coils, N_meas]
        
        inds = np.argsort(np.concatenate((twix_obj.image.memPos, twix_obj.phasecor.memPos)))
        
        kdata = np.concatenate((kdata, kdata_phasecor), axis=-1)[:,:,inds]
    
    # Reorder to: [N_coils, N_meas, N_adc]
    kdata = kdata.transpose(1,2,0)
    
    return kdata

# Sum of squares coil combination, assumes x is a matrix of size [N_coils, ...]
def combine_coils(x):
    return np.sqrt((abs(x)**2).sum(axis=0))

def ifft_2d(x):
    return ifft_nd(x, (-2,-1))

def ifft_nd(x, axes):
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x, axes=axes), axes=axes), axes=axes)

def sort_data_labels(kdata, seq, shape=None):
    n_coils = kdata.shape[0]
    adc_samples = kdata.shape[2]

    supported_labels = ['SEG', 'SET', 'REP', 'PHS', 'ECO', 'SLC', 'LIN', 'AVG']

    # Get label evolutions from sequence
    labels = seq.evaluate_labels(evolution='adc')

    # Reverse lines
    if 'REV' in labels:
        rev = labels['REV'] != 0
        kdata[:, rev, :] = kdata[:, rev, ::-1]

    # Find out which labels are used in the sequence and calculate kspace shape
    index = []
    label_shape = []
    used_labels = []
    for lab in supported_labels:
        if lab in labels:
            index.append(labels[lab])
            label_shape.append(labels[lab].max() + 1)
            used_labels.append(lab)
    
    label_shape += [adc_samples]
    
    if shape is None:
        shape = label_shape
        print(f'Automatically detected matrix size: {used_labels + ["ADC"]} {shape}')
    elif len(shape) != len(label_shape):
        raise ValueError('Provided shape does not have the same number of dimensions as the number of labels')

    # Assigned measured data to kspace
    kspace_matrix = np.zeros((n_coils,) + tuple(shape), dtype=np.complex64)
    kspace_matrix[(slice(None),) + tuple(index) + (slice(None),)] = kdata

    # Do averaging
    if 'AVG' in labels:
        kspace_matrix = kspace_matrix.mean(axis=-2)

    return kspace_matrix


def recon_cartesian_2d(kdata, seq, shape=None, use_labels=None):
    if use_labels is None:
        # Detect if sequence used labels
        use_labels = seq.label_set_library.data != {} or seq.label_inc_library.data != {}
        if use_labels:
            print('Detected labels in the sequence!')
        else:
            print('Did not detect labels in the sequence, using kspace calculation for sorting!')
    
    if use_labels:
        kspace = sort_data_labels(kdata, seq, shape=shape)
    else:
        kspace = sort_data_implicit(kdata, seq, shape=shape)
    im = ifft_2d(kspace)
    if im.shape[0] > 1:
        sos = combine_coils(im)
    else:
        sos = im[0]
    
    return sos

def main():
    ### READ DATA ###
    parser = ArgumentParser(description='Reconstruct raw data from pulseq mese')
    parser.add_argument('-p', '--path', type=str, default= os.getcwd(), help='path to the dataset, default: current working directory')
    parser.add_argument('-a', '--anonymize', const='anon', metavar='pseudo_name', dest='anonymize', type=str, nargs = '?', help='Use the pseudo_name (default: anon) as patient name')
    args = parser.parse_args()

    path = args.path
    ANON_NAME = args.anonymize # not yet used, have to figure out how to best do this with multiple .dat files
    
    # get file names for Pulseq data
    all_files=os.listdir(path)
    
    seq_filename = [x for x in all_files if x.endswith('.seq')]
    if seq_filename == []:
        raise ValueError('No .seq file uploaded')
    seq_filename = seq_filename[0]
    
    filenames = [x for x in all_files if x.endswith('.dat')]
    if filenames == []:
        raise ValueError('No .dat file uploaded')
    
    json_filenames = [x for x in all_files if x.endswith('.json')]
    json_filename = []
    for i in range(0,len(json_filenames)):
        if re.search('(.+?)[.]seq', seq_filename).group(1) == re.search('(.+?)[.]json', json_filenames[i]).group(1):
            json_filename = json_filenames[i]
    
    if json_filename == []:
        raise ValueError('No matching .json file found')
    
    for i in range(0, len(filenames)):
        filename = filenames[i]
        # Load data from twix file (shape = [N_coils, N_meas, N_adc])
        np.seterr(divide='ignore', invalid='ignore')
        kdata = read_raw_data(os.path.join(path, filename))
    
        # Load associated sequence file
        seq = Sequence()
        seq.read(os.path.join(path, seq_filename), detect_rf_use=True)
    
        # Load json file
        with open(os.path.join(path, json_filename)) as json_file:
            jsonfile = json.load(json_file)
            json_file.close()
    
        # get header
        twix = tx.read_twix(os.path.join(path, filename))
        hdr = twix[-1]['hdr']
    
        # get header and json info
        pos_pat = twix[-1]['hdr']['Config']['PatientPosition']
        freq = twix[-1]['hdr']['Dicom']['lFrequency']
        field_strength = twix[-1]['hdr']['Dicom']['flMagneticFieldStrength']
        res = jsonfile['Resolution']
        slice_thickness = jsonfile['SliceThickness']
    
        # reconstruct data
        rec = recon_cartesian_2d(kdata, seq, use_labels=True) # combines coils internally
        # reorder pulseq data as [x,y,z,echoes], scale it correctly and make nifti file 
        rec_reordered = np.transpose(rec, [3,2,0,1])
        rec_reordered -= np.min(rec_reordered)
        rec_reordered *= 4096/np.max(rec_reordered)
        dims = np.shape(rec)
        dims_re = np.shape(rec_reordered)
        
        new_dir = pathlib.Path(path, 'mr-anat')
        new_dir.mkdir(parents=True, exist_ok=True)
    
        new_filename = re.search('(.+?)[.]dat', filename).group(1) + '_mese'
        nii_image = nib.Nifti1Image(rec_reordered, affine=np.eye(4)*[res[0], res[1], slice_thickness, 1])
        nib.save(nii_image, os.path.join(path, new_dir.name, new_filename + '.nii.gz'))
    
        # write json file for this data
        additional_entries = {'PatientPosition': pos_pat,
                             'MagneticFieldStrength': field_strength,
                             'ImagingFrequency': freq}
        jsonfile.update(additional_entries)
    
        with open(os.path.join(path, new_dir.name, new_filename + '.json'), 'w') as f:
                json.dump(jsonfile, f, indent=2)
            
            
if __name__== '__main__':
    main()
