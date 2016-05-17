*The MS3 folder refers to the folder downloaded from https://github.com/JessCG/MS3.git
It may be useful to rename the folder after download and move it to a convenient 
location on your computer.

The MS3 folder contains all necessary functions and folder structure to run the script
MergeProcessor.m (MS3 > matlab). The script was written by John Newgard 
(Dalhousie University) to merge size distributions obtained with a Multisizer 3, using 
2 aperture tubes (30 and 200 microns) or 3 aperture tubes (30, 200 and 400 microns). 
When dilution sheets are included, the script merges using PPM values. When no dilution 
sheets are included, the script merges using % volume. More information can be found in 
the script itself.

In order to run the script, users must enter the location of the MS3 folder on their
computer in the file mergedir.m (MS3 > matlab).

Users should also modify the MS3_flowrates.txt file found in the folder MS3 by 
entering flow rates specific to their instrument.

Questions regarding the download should be directed to Jessica (jgarwood@ucsd.edu).

---

Folder descriptions:

-Dilutions: dilution sheets (follow sample format) corresponding to files to be merged 
should be put in this directory.

-ToMergeProcessor: MS3 data output (in .CSV format) should be put in this folder. All files 
within the folder will be merged so it is important to remove previously merged files. 
Sample output formats are included with this download.

-MergedData: merged distributions will be saved in this location.

-matlab: MergeProcessor.m script and functions necessary to run this script
