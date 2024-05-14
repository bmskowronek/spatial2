# spatial2

# recommended guide through these files is in raport.pdf

Note: the mentioned big plots can be found here: https://drive.google.com/file/d/12j5eDyuZek_BXfgjw7pSenfhk3HRZUs3/view?usp=sharing
However, mind that the code plot_app7_raw.py should be able to generate them anyway, albeit a bit slow.
Note: this version may have some small bugs

All graphing files used a modified, reduced version of the data, transformed by the file "datachanger.py". 

Run this datachanger.py file first 
WARNING: it will work only on dataset that has files that are missing phenotypes removed! 
warning2: using it should engage about 632MB of additional disk space

with a correct path set in the script to where your if_data folder is placed, to continue using the other files with plots, graphs etc. This is due to unintelligible performance issues caused by attempting to import dataframes and other content from one python file to another while developing this. Most plotting files generate results in external html files. I considered developing this as a streamlit application which might have been way easier to use, however various performance issues made me tentatively reject it and it remained rejected until this version.
