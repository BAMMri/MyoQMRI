# MyoQMRI

This project is an open-source effort to put together tools for quantitative MRI 
of the muscles. 
Current features include:

* multi echo spin echo sequence, written in PyPulseq
* fast water T2 mapping from multi echo spin echo images
* fat water separation from 2 echo gradient echo data

## Setup

### Virtual environment

We recomment creating a virtual environment.  
Do this with

```
python -m venv /path/to/new/virtual/environment
```

If you are on Windows, this will look something like this:

```
python -m venv C:\path\to\new\virtual\environment
```

Depending on your distribution, you may need to use python3 instead of python as
command.  
If there are any problems concerning the creation of the virtual environment or
its use, refer to:
https://docs.python.org/3/library/venv.html

To activate the virtual environment, run:

```
source /path/to/new/virtual/environment/bin/activate
```

On Windows:

```
C:\path\to\new\virtual\environment\Scripts\activate.bat
```

To deactivate the virtual environment, run:

```
deactivate
```

### Installation

#### Install via PyPI

MyoQMRI is on PyPI, you can easily install it via 

```
pip install myoqmri
```
if you have cuda, also install the optional dependencies

```
pip install myoqmri[cuda]
```

#### Download and install via GitHub

Alternatively, you can install MyoQMRI directly from GitHub.
The git repository contains the compiled pulseq multi echo spin echo sequence, 
all analysis tools and example data.  
Download it via

```
git clone https://github.com/bammri/MyoQMRI
```

or manually via the webpage.  
Install it with

```
pip install -e /path/to/folder/MyoQMRI
```

If you want to use cuda for GPU support, install the optional dependencies via

```
pip install -e /path/to/folder/MyoQMRI[cuda]
```

On Windows:

```
pip install -e C:\path\to\folder\MyoQMRI
```


## Multi-Echo Spin Echo Sequence

For instructions on how to install the pulseq interpreter and bring the sequence 
to the scanner, see 
[https://github.com/pulseq/tutorials](https://github.com/pulseq/tutorials).  

The sequence was developed for imaging of the thighs. In principle, however, 
it is possible to easily adapt it for other body regions.  
For this,

* Copy the `write_pulseq_mese_leg.py` file
* Rename the file and the `seq_filename` variable within
* Edit the `fov_x` and `fov_y` variables (and the `Nx, Ny` accordingly to keep 
the resolution the same)
* Execute the code in the terminal with: `python write_pulseq_your_seqname.py`

The `.seq` and `.json` file are written into the current working directory.

### Reconstruction

Data from the `mese_leg.seq` is reconstructed offline.

* Transfer the raw data from the scanner
* rename the `.dat` file to the patient token of your choice
* Copy `mese_leg.json` and `mese_leg.seq` into the same folder
* Run `raw2nii_pulseq_mese_leg -p path/to/folder/containing/the/data`

The code outputs a `.nii.gz` and an according `.json` file into an `/mr-anat` folder.

<a name="prepare-dicom-data">

## Prepare DICOM data
</a>

Dicom files need to be converted to the ormir-mids format.  
For this, we use the ormir-mids converter (more information on 
[https://github.com/ormir-mids/ormir-mids](https://github.com/ormir-mids/ormir-mids)).
If you are working with the enhanced dicom format and are converting Dixon data, 
see additional [step below](#additional-step-for-enhanced-dicom).

The conversion is done with:

```
dcm2omids -r -a PATIENT_NAME /path/to/input/folder /path/to/output/folder 
```

Note: Output folder must exist.

<a name="additional-step-for-enhanced-dicom">

#### Additional step for enhanced dicom
</a>

If you are working with the old dicom format, you can ignore this step.

As of now, for enhanced Dixon dicoms, we need a series_config.json file for the 
converter to work properly.
For this, each magnitude and phase 4D image needs to be in its own folder (see 
example data folder structure).
Copy the series_config.json file from the example data and edit it to 
contain the according series number of your data. Afterwards, the conversion can 
be executed like above.

## T2 Mapping From Multi-Echo Spin Echo

The main function is used to obtain **water T2** and **fat fraction** from
 multi-echo spin echo images as described in 
 [Marty et al. ](https://doi.org/10.1002/nbm.3459)

This implementation uses the GPU for the generation of a dictionary of 
signals through [Extended Phase Graph simulation](https://doi.org/10.1002/jmri.24619).

Ideally, you should work with the ormir-mids format. See Section 
[Prepare DICOM data](prepare-dicom-data) above for instructions. However, as of now, it is also possible 
to use the code with dicom images.

In case of dicom data, the dataset must be a directory containing 2D DICOM 
images from a multiecho spin echo acquisition, ordered as 
(slice1echo1 - slice2echo1 - ... - slice1echo2 - slice2echo2 - ...).

### Usage

   usage: waterT2 [-h] [--fit-type T] [--bids] [--path-is-nifti] [--fat-t2 T2] 
                  [--noise-level N] [--nthreads T] [--plot-level L] 
                  [--t2-limits min max] [--b1-limits min max] [--use-gpu] 
                  [--ff-map dir] [--register-ff] [--etl-limit N] 
                  [--out-suffix ext] [--slice-range start end] 
                  [--refocusing-width factor] [--exc-profile path] 
                  [--ref-profile path] path
    
    Fit a multiecho dataset
    
    positional arguments:
      path                  path to the dataset
    
    optional arguments:
      --bids                use muscle-BIDS format for input/output
      --path-is-nifti       set if path is pointing directly to a nifti file
      -h, --help            show this help message and exit
      --fit-type T, -y T    type of fitting: T=0: EPG, T=1: Single exponential, T=2: Double exponential (default: 0)
      --fat-t2 T2, -f T2    fat T2 (default: 151)
      --noise-level N, -n N
                            noise level for thresholding (default: 300)
      --nthreads T, -t T    number of threads to be used for fitting (default: 12)
      --plot-level L, -p L  do a live plot of the fitting (L=0: no plot, L=1: show the images, L=2: show images and signals)
      --t2-limits min max   set the limits for t2 calculation (default: 20-80)
      --b1-limits min max   set the limits for b1 calculation (default: 0.5-1.2)
      --use-gpu, -g         use GPU for fitting
      --ff-map dir, -m dir  load a fat fraction map
      --register-ff, -r     register the fat fraction dataset
      --etl-limit N, -e N   reduce the echo train length
      --out-suffix ext, -s ext
                            add a suffix to the output map directories
      --slice-range start end, -l start end
                            Restrict the fitting to a subset of slices
      --refocusing-width factor, -w factor
                            Slice width of the refocusing pulse with respect to the excitation (default 1.2) (Siemens standard)
      --exc-profile path    Path to the excitation slice profile file
      --ref-profile path    Path to the refocusing slice profile file

#### Example

```
waterT2 --bids --path-is-nifti -n100 PATNAME_mese.nii.gz
```

### Slice profile

An accurate slice profile is crucial to obtain accurate results. By default,
a hanning-windowed sinc pulse is used, and the refocusing pulse has a 1.2x the
slice width of the excitation pulses. This reflects the parameters of the
standard Siemens Spin Echo sequence.

External slice profile files can be provided. They are text files with angle
values (in degrees) across half the profile (i.e. starting with ~90 or ~180 and
decreasing), with one value per line.

Either none or both slice profiles must be given, and both slice arrays must
contain the same number of samples.

Example:

    90
    89.8014
    89.2106
    88.2423
    86.9198
    85.2734
    83.3385
    81.1537
    78.7585
    76.1921
    ...


## Fat Water Separation
The 2 echo fat water separation uses an algorithm by the BMRR group in 
Munich, which utilizes hirarchical multi-resolution graph-cuts (for more 
information see 
[https://github.com/BMRRgroup/2echo-WaterFat-hmrGC](https://github.com/BMRRgroup/2echo-WaterFat-hmrGC)).  

### Download and Install Package

Download the code via

```
git clone https://github.com/BMRRgroup/2echo-WaterFat-hmrGC
```

or manually via the webpage.  
Install it with

```
pip install -e /path/to/folder/2echo-WaterFat-hmrGC
```

### Usage

    usage: fatwater [-h] [-p PATH] [-e ECHONUMS ECHONUMS] [-c FATSHIFT [FATSHIFT ...]] [-a RELAMPS [RELAMPS ...]] [-ph]

    Compute fat/water images from dual echo data

    options:
      -h, --help            show this help message and exit
      -p PATH, --path PATH  path to the folder where the data is located, 
                            default: current working directory
      -e ECHONUMS ECHONUMS, --echonums ECHONUMS ECHONUMS
                            if data contains more than two echoes, provide which
                            echoes to use, default: first two echoes
      -c FATSHIFT [FATSHIFT ...], --fatshift FATSHIFT [FATSHIFT ...]
                            Chemical shift of fat peak(s) in ppm, default: 3.5
      -a RELAMPS [RELAMPS ...], --relamps RELAMPS [RELAMPS ...]
                            Relative amplitude of fat peaks, default: 1
      -ph                   if given, the code additionally outputs the phase of
                            the computed fat and water images

#### Example

```
fatwater -p path/to/folder
```

## Git repositories of the constituents

### Pulseq/PyPulseq

https://pulseq.github.io/index.html  

### Ormir-mids

https://github.com/ormir-mids

### MyoQMRI

https://github.com/bammri/MyoQMRI

### Fat water separation

https://github.com/BMRRgroup/2echo-WaterFat-hmrGC


## Citing

### PyPulseq

Ravi KS, Geethanath S, Vaughan JT. PyPulseq: A python package for MRI pulse 
sequence design. Journal of Open Source Software, 2019;4:1725.  

### Pulseq mese sequence

Schäper J, Santini F, Weidensteiner C. *A standardized open-source MESE sequence 
implemented in PyPulseq for reproducible water T2 quantification in skeletal 
muscle*, Proceedings of the International Society for Magnetic Resonance in 
Medicine (ISMRM) 2025  

### Water T2 quantification

Santini F, Deligianni X, Paoletti M, et al., *Fast open-source toolkit for water 
T2 mapping in the presence of fat from multi-echo spin-echo acquisitions for 
muscle MRI*, Frontiers in Neurology 2021, https://doi.org/10.3389/fneur.2021.630387.  

### Fat water separation

Stelter J et al. *Hierarchical multi-resolution graph-cuts for water-fat-silicone 
separation in breast MRI*, IEEE Transactions on Medical Imaging, 2022, 
DOI: 10.1109/TMI.2022.3180302, https://ieeexplore.ieee.org/document/9788478  

Eggers H, Brendel B, Duijndam A, Herigault G. *Dual-echo Dixon imaging with 
flexible choice of echo times*, Magnetic Resonance in Medicine, 2010, 
DOI: 10.1002/mrm.22578, https://doi.org/10.1002/mrm.22578  

# Authors and Acknowledgment

### Main contributors

* Jessica Schäper
* Francesco Santini
* Claudia Weidensteiner

EPG code was made possible by Matthias Weigel.

This project was supported by

* [FSRMM](https://www.fsrmm.ch/)
* [SNF](http://www.snf.ch/) (grant number 320030_172876)







