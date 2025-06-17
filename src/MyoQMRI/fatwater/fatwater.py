"""
    This code computes fat and water images from 2 echo GRE data, using code
    which was developed by Jonathan Stelter et al.
    For details see:
    https://github.com/BMRRgroup/2echo-WaterFat-hmrGC for instructions
    Optional input arguments:
    -p --path       path to folder which contains data (default: current directory)
    -e --echonums   if data contains more than two echoes, provide which echoes
                    to use (default: first two echoes)
    -c --fatshift   Chemical shift of fat peak(s) in ppm (default: 3.5)
    -a --relamps    (Relative amplitude of fat peaks (default: 1), required if
                    more than one argument was given for the fatshift
    -ph             if given, the code additionally outputs the phase of the 
                    computed fat and water images
"""

import os
import re
import pathlib
import numpy as np
from argparse import ArgumentParser
import nibabel as nib
import json

try:
    from hmrGC_dualEcho.dual_echo import DualEcho
except ImportError:
    raise RuntimeError('hmrGC_dualEcho is not available. Visit https://github.com/BMRRgroup/2echo-WaterFat-hmrGC for instructions')


def main():
    #### READ ARGUMENTS ###
    parser = ArgumentParser(description='Compute fat/water images from dual echo data')
    parser.add_argument('-p', '--path', type=str, default= os.getcwd(), help='path to the folder where the data is located, default: current working directory')
    parser.add_argument('-e', '--echonums', nargs=2, type=int, default=[1,2], help='if data contains more than two echoes, provide which echoes to use, default: first two echoes')
    parser.add_argument('-c', '--fatshift', nargs='+', type=float, default=[3.5], help='Chemical shift of fat peak(s) in ppm, default: 3.5')
    parser.add_argument('-a', '--relamps', nargs='+', type=float, default=[1], help='Relative amplitude of fat peaks, default: 1')
    parser.add_argument('-ph', action='store_true', help='if given, the code additionally outputs the phase of the computed fat and water images')
    args = parser.parse_args()

    path = args.path
    eno1 = args.echonums[0]
    eno2 = args.echonums[1]
    fatshift = args.fatshift
    relamps = args.relamps

    if len(fatshift) != len(relamps):
        raise ValueError('Number of fat peaks and relative amplitudes must match')

    # get file names
    all_files = os.listdir(path)
    magnnames = [x for x in all_files if x.endswith('megre.nii.gz')]
    phasenames = [x for x in all_files if x.endswith('megre_ph.nii.gz')]
    jsonfilename = [x for x in all_files if x.endswith('megre.json')]
    if magnnames == []:
        raise ValueError('No magnitude data found')
    if phasenames == []:
        raise ValueError('No phase data found')
    if jsonfilename == []:
        raise ValueError('No json file found')

    magnnames.sort()
    phasenames.sort()
    jsonfilename.sort()

    ### PREPARE DATA, EXECUTE FAT WATER COMPUTATION, SAVE AS NIFTI ###
    for i in range(0, len(magnnames)):
    # load data as list, convert to np array
        magndata = nib.load(os.path.join(path, magnnames[i]))
        print(magnnames[i] + ' loaded')
        affinematrix = magndata.affine
        pixeldim = magndata.header['pixdim']
        pixeldim = pixeldim[1:4]
        magndata = magndata.get_fdata()
        phasedata = nib.load(os.path.join(path,phasenames[i]))
        print(phasenames[i] + ' loaded')
        phasedata = phasedata.get_fdata()
        phasedata = (phasedata - 2048)/4096 * 2 * np.pi
    
    # combine data to complex array
        data_complex = magndata * np.exp(1j*phasedata)
    
    # load json file and get relevant parameters
        with open(os.path.join(path, jsonfilename[i])) as json_file:
            jsonfile = json.load(json_file)
            json_file.close()
    
        field_strength = jsonfile['MagneticFieldStrength']
        center_freq = jsonfile['ImagingFrequency']*(10**6)
        TEs = jsonfile['EchoTime'] 
    
    # calculate mask from magnitude
        msk = abs(magndata[:,:,:,1])
        msk = msk/msk.max()
        msk[msk < 0.02] = 0
        msk[msk > 0] = 1
    
    ## get 2 specific echos
        TE1 = TEs[eno1-1]
        TE2 = TEs[eno2-1]
        TEs = [TE1, TE2]
        data_complex1 = data_complex[:,:,:,eno1-1]
        data_complex2 = data_complex[:,:,:,eno2-1]
        data_complex = np.stack((data_complex1, data_complex2), axis=3)
    
    # Input arrays and parameters
        signal = data_complex   # complex array with dim (nx, ny, nz, nte)
        mask = msk   # boolean array with dim (nx, ny, nz)
        params = {}
        params['TE_s'] = np.asarray(TEs)*(10**(-3))   # float array with dim (nte)
        params['centerFreq_Hz'] = center_freq   # float (in Hz, not MHz)
        params['fieldStrength_T'] = field_strength   # float
        params['voxelSize_mm'] = np.asarray(pixeldim)   # recon voxel size with dim (3)
        params['FatModel'] = {}
        params['FatModel']['freqs_ppm'] = np.asarray(fatshift)   # chemical shift difference between fat and water peak, float array with dim (nfatpeaks)
        params['FatModel']['relAmps'] = np.asarray(relamps)   # relative amplitudes for each fat peak, float array with dim (nfatpeaks)
    
    # Initialize DualEcho object
        g = DualEcho(signal, mask, params)
    
    # Perform graph-cut method
        g.perform()   # methods with different parameters can be defined using the dual_echo.json file
    
    # access separation results
        watermagn = abs(g.images['water'])
        fatmagn = abs(g.images['fat'])
        waterphase = np.angle(g.images['water'])
        fatphase = np.angle(g.images['fat'])
    
    # export as nifti
        try:
            patname = re.search('(.+?)_megre', magnnames[i]).group(1)
        except AttributeError:
            patname = '' 
    
        new_dir = pathlib.Path(path, 'mr-quant')
        new_dir.mkdir(parents=True, exist_ok=True)
        filename_water = patname + '_W.nii.gz'
        filename_fat = patname + '_F.nii.gz'
    
        nii_image = nib.Nifti1Image(watermagn, affine=affinematrix)
        nib.save(nii_image, os.path.join(path, new_dir.name, filename_water))
        nii_image = nib.Nifti1Image(fatmagn, affine=affinematrix)
        nib.save(nii_image, os.path.join(path, new_dir.name, filename_fat))
        if args.ph == True:
            filename_water_phase = patname + '_W_phase.nii.gz'
            filename_fat_phase = patname + '_F_phase.nii.gz'
            nii_image = nib.Nifti1Image(waterphase, affine=affinematrix)
            nib.save(nii_image, os.path.join(path, new_dir.name, filename_water_phase))
            nii_image = nib.Nifti1Image(fatphase, affine=affinematrix)
            nib.save(nii_image, os.path.join(path, new_dir.name, filename_fat_phase))

        # write json file for this data
        additional_entries_water = {'PulseSequenceType': 'Water Map'}
        additional_entries_fat = {'PulseSequenceType': 'Fat Map'}
        jsonfile.update(additional_entries_water)
        jsonfile_water = jsonfile
        jsonfile.update(additional_entries_fat)
        jsonfile_fat = jsonfile
        with open(os.path.join(path, new_dir.name, patname + '_W.json'), 'w') as f:
                json.dump(jsonfile_water, f, indent=2)
        with open(os.path.join(path, new_dir.name, patname + '_F.json'), 'w') as f:
                json.dump(jsonfile_fat, f, indent=2)
           
if __name__ == '__main__':
    main()
