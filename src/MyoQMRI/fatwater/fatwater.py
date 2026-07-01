"""
    This script computes fat and water images from 2-echo GRE data, using code
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
from pathlib import Path
from collections import defaultdict
import numpy as np
from argparse import ArgumentParser
import nibabel as nib
import json
import copy

try:
    from hmrGC_dualEcho.dual_echo import DualEcho
except ImportError:
    raise RuntimeError('hmrGC_dualEcho is not available. Visit https://github.com/BMRRgroup/2echo-WaterFat-hmrGC for instructions')

##################### Functions for reading data ###############################
def find_megre_files(folder: str | Path) -> dict:
    """
    Scan a folder for MEGRE files and group them by dataset prefix.

    Returns a dict keyed by prefix (e.g. "xxx", "yyy"), each containing
    'magnitude' and 'phase' sub-dicts with 'json' and 'nii' paths.
    """
    folder = Path(folder)

    # prefix -> {"magnitude": {"json":..., "nii":...}, "phase": {...}}
    datasets = defaultdict(lambda: {
        "magnitude": {"json": None, "nii": None},
        "phase":     {"json": None, "nii": None},
    })

    phase_pattern = re.compile(r"^(.*?)_part-phase_MEGRE(\..+)$")
    mag_pattern   = re.compile(r"^(.*?)(?:_part-mag)?_MEGRE(\..+)$")

    for file in folder.iterdir():
        name = file.name

        if name.endswith(".json"):
            ext_key = "json"
        elif name.endswith(".nii.gz") or name.endswith(".nii"):
            ext_key = "nii"
        else:
            continue

        m_phase = phase_pattern.match(name)
        if m_phase:
            prefix = m_phase.group(1)
            datasets[prefix]["phase"][ext_key] = file
            continue

        m_mag = mag_pattern.match(name)
        if m_mag:
            prefix = m_mag.group(1)
            datasets[prefix]["magnitude"][ext_key] = file

    return dict(datasets)


def load_megre_data(folder: str | Path) -> dict:
    """
    Load MEGRE JSON metadata + nii paths for every dataset found in the folder,
    keyed by dataset prefix.
    """
    all_files = find_megre_files(folder)
    all_data = {}

    for prefix, files in all_files.items():
        data = {}
        for kind in ("magnitude", "phase"):
            json_path = files[kind]["json"]
            nii_path  = files[kind]["nii"]

            data[kind] = {
                "metadata": json.loads(json_path.read_text()) if json_path else None,
                "nii_path": nii_path,
            }

            status = []
            status.append(f"JSON={json_path.name}" if json_path else "JSON=missing")
            status.append(f"NII={nii_path.name}" if nii_path else "NII=missing")
            print(f"[{prefix} | {kind:>9}] {' | '.join(status)}")

        all_data[prefix] = data

    return all_data

################################################################################

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
    ph = args.ph

    if len(fatshift) != len(relamps):
        raise ValueError('Number of fat peaks and relative amplitudes must match')

    #load data
    data = load_megre_data(path)

    ### PREPARE DATA, EXECUTE FAT WATER COMPUTATION, SAVE AS NIFTI ###
    keys = list(data)
    for i in range(0, len(keys)):
    # load data as list, convert to np array
        magndata = nib.load(data[keys[i]]["magnitude"]["nii_path"])
        print(str(data[keys[i]]["magnitude"]["nii_path"]) + ' loaded')
        affinematrix = magndata.affine
        pixeldim = magndata.header['pixdim']
        pixeldim = pixeldim[1:4]
        magndata = magndata.get_fdata()
        phasedata = nib.load(data[keys[i]]["phase"]["nii_path"])
        print(str(data[keys[i]]["phase"]["nii_path"]) + ' loaded')
        phasedata = phasedata.get_fdata()
        phasedata = (phasedata - 2048)/4096 * 2 * np.pi
    # combine data to complex array
        data_complex = magndata * np.exp(1j*phasedata)
    
    # access relevant meta data
        field_strength = data[keys[i]]["magnitude"]['metadata']['MagneticFieldStrength']
        center_freq = data[keys[i]]["magnitude"]['metadata']['ImagingFrequency']*(10**6)
        TEs = data[keys[i]]["magnitude"]['metadata']['EchoTime'] 
    
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
        patname = keys[i]
    
        new_dir = Path(path, 'mr-quant')
        new_dir.mkdir(parents=True, exist_ok=True)
        if ph == True:
            filename_water = patname + '_part-mag_WATER'
            filename_fat = patname + '_part-mag_FAT'
            filename_water_phase = patname + '_part-phase_WATER'
            filename_fat_phase = patname + '_part-phase_FAT'
            nii_image = nib.Nifti1Image(waterphase, affine=affinematrix)
            nib.save(nii_image, os.path.join(path, new_dir.name, filename_water_phase + '.nii.gz'))
            nii_image = nib.Nifti1Image(fatphase, affine=affinematrix)
            nib.save(nii_image, os.path.join(path, new_dir.name, filename_fat_phase + '.nii.gz'))
        else:
            filename_water = patname + '_WATER'
            filename_fat = patname + '_FAT'
    
        nii_image = nib.Nifti1Image(watermagn, affine=affinematrix)
        nib.save(nii_image, os.path.join(path, new_dir.name, filename_water + '.nii.gz'))
        nii_image = nib.Nifti1Image(fatmagn, affine=affinematrix)
        nib.save(nii_image, os.path.join(path, new_dir.name, filename_fat + '.nii.gz'))
    
        # write json file for this data
        metadata_water = copy.deepcopy(data[keys[i]]["magnitude"]['metadata'])
        metadata_fat = copy.deepcopy(data[keys[i]]["magnitude"]['metadata'])
        metadata_water['PulseSequenceType'] = 'Water Map'
        metadata_fat['PulseSequenceType'] = 'Fat Map'
        with open(os.path.join(path, new_dir.name, filename_water + '.json'), 'w') as f:
                json.dump(metadata_water, f, indent=2)
        with open(os.path.join(path, new_dir.name, filename_fat + '.json'), 'w') as f:
                json.dump(metadata_fat, f, indent=2)
        if ph == True:
            with open(os.path.join(path, new_dir.name, filename_water_phase + '.json'), 'w') as f:
                    json.dump(metadata_water, f, indent=2)
            with open(os.path.join(path, new_dir.name, filename_fat_phase + '.json'), 'w') as f:
                    json.dump(metadata_fat, f, indent=2)

if __name__ == '__main__':
    main()