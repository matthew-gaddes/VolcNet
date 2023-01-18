# VolcNet
A database of labelled Sentinel-1 data featuring examples of volcanic unrest.  




# Version 3.X 
The first version of this is due in early 2023 and will feature a more accisible file format and a significant increase in the number of volcanoes.  


# Version 2.X 
Features the data that can be used to create ~500,000 labelled interferograms.  The updates are detailed in version three of "Simultaneous classification and location of volcanicdeformation in SAR interferograms using deep learningand the VolcNet database" (https://eartharxiv.org/repository/view/1969/), but to summarise we swtich to using time series with labels of the duration and magnitude of deformation, and then automatically create labelled interferograms between any two acquisitions in the time series.  

To create your own labelled data, have a a look at the example in bin/02_volcnet_labelled_ifg_constructor.py

An overview of how all possible interferograms between all acquisitions can be made and labelled for Sierra Negra.  

![figure_7_volcnet_sierra_negra](https://user-images.githubusercontent.com/10498635/213170308-f43892c3-e411-4df0-a651-d239d55e9e8a.png)

A breakdown of the type and label of data:
![figure_6_volcnet_summary](https://user-images.githubusercontent.com/10498635/213170457-5dd87f87-e332-4d7d-b0ae-ad2d941cf94e.png)





# Version 1.0 
Features a set of ~250 labelled unwrapped interferograms that contain labels of both the type of deformation (including examples of no deformation), and the location of deformation within the interferograms.  These have been used to train a CNN, VUDL-21, which was detailed in the first two vesion of "Simultaneous classification and location of volcanicdeformation in SAR interferograms using deep learningand the VolcNet database" (https://eartharxiv.org/repository/view/1969/).  The code for the CNN is also on [Github](https://github.com/matthew-gaddes/VUDLNet_21).  

The figures below show examples of the data and labels.  For simplicity, the database requires only pickle to open the files, and numpy to handle arrays, and interferograms are stored using masked arrays in order for incoherent pixels to be masked.  As many interferograms are from time series, these are collated into a single file, and each file may therefore be of different lengths.  V1.1 plans to include an index function to find interferograms given a volcanoes name.  

![02_Real_data](https://user-images.githubusercontent.com/10498635/104197662-37177100-541d-11eb-9cc5-749ec6ea5caa.png)
![02_Real_data2](https://user-images.githubusercontent.com/10498635/104197668-37b00780-541d-11eb-9f36-2054fb91cb77.png)
![02_Real_data3](https://user-images.githubusercontent.com/10498635/104197672-38489e00-541d-11eb-8843-a966a400d9ca.png)

