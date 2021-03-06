# THINGSanalysis

Data and code to analyse the luminance and contrast of images in the [THINGS database](https://osf.io/jum2f/) as in my 
recent work published in [Perception](https://journals.sagepub.com/eprint/9CPFR9QQ5JJPYDFEUABR/full) also available via  [bioRXiv](https://www.biorxiv.org/content/10.1101/2021.07.08.451706v2).

**Harrison, W. J.** 2022. “Luminance and Contrast of Images in the THINGS Database:” Perception 51(4):244–62. doi: 10.1177/03010066221083397.

![Figure 1](./output/allIms.png?raw=true)

## Table of Contents
* [Installation and requirements](#installation-and-requirements)
* [Instructions for the code](#instructions-for-the-code)
* [Citation](#citation)
* [Contact](#contact)

## Installation and requirements

- Clone this repository:
    
    
    git clone https://github.com/willjharrison/THINGSanalysis.git
    
    
- Download [THINGS dataset](https://osf.io/jum2f/). 
    - Extract the dataset (zip file) to this directory (path_to/THINGSanalysis);
    - Rename the folder to 'THINGS'
    - Extract the zip files in 'THINGS/Main/' to 'THINGS/Main/images/'. Password is available here [THINGS dataset](https://osf.io/jum2f/)

## Instructions for the code
All analysis code are Matlab files, and all images and data generated by the code are also given here. Images and data files are sent to the "output" folder, so you can use that if you don't want to run the analyses yourself (some of them take a bit of time).

THINGS_imageStatistics_v02.m - analyses statistics <- this needs to be run before anything else will work

classification_v01.m - classifies images in pairs according to leave one out analysis

classification_allCats_v01.m - classifies images into one of 1854 concepts, as opposed to pairs

classification_allStats_v01.m - classifies images in pairs using spatial frequency and oriented energy in addition to mean luminance and RMS contrast

THINGS_imageStatistics_v01_normedIms_v01.m - repeats the pairwise classification, but for normalised images

Data figures in the manuscript were created with DataGraph.

**If you're worried about these confounds in your own experiment, you can use/edit the following function to normalise the contrast and luminance while retaining the colour of your stimui.**

normColIms.m - shows an example of how to normalise the luminance and RMS contrast of coloured images

## Citation
    
    @article{Harrison2022,
        title = {Luminance and {Contrast} of {Images} in the {THINGS} {Database}:},
        volume = {51},
        copyright = {© The Author(s) 2022},
        url = {https://journals.sagepub.com/eprint/9CPFR9QQ5JJPYDFEUABR/full},
        doi = {10.1177/03010066221083397},
        number = {4},
        journal = {Perception},
        author = {Harrison, William J.},
        year = {2022},
        note = {Publisher: SAGE PublicationsSage UK: London, England},
        pages = {244--262},
        
    }


## Contact

Will Harrison <[willjharri@gmail.com](willjharri@gmail.com)>
