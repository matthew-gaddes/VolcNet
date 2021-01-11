# VolcNet
A database of labelled Sentinel-1 data featuring examples of volcanic unrest.  Version 1.0 features a set of ~250 labelled unwrapped interferograms that contain labels of both the type of deformation (including examples of no deformation), and the location of deformation within the interferograms.  These have been used to train a convolutional neural network, VUDL-21, which can also be found on [Github](https://github.com/matthew-gaddes/VUDLNet_21).  

The figures below show examples of the data and labels.  For simplicity, the database requires only pickle to open the files, and numpy to handle arrays, and interferograms are stored using masked arrays in order for incoherent pixels to be masked.  As many interferograms are from time series, these are collated into a single file, and each file may therefore be of different lengths.  V1.1 plans to include an index function to find interferograms given a volcanoes name.  

![02_Real_data](https://user-images.githubusercontent.com/10498635/104197662-37177100-541d-11eb-9cc5-749ec6ea5caa.png)
![02_Real_data2](https://user-images.githubusercontent.com/10498635/104197668-37b00780-541d-11eb-9f36-2054fb91cb77.png)
![02_Real_data3](https://user-images.githubusercontent.com/10498635/104197672-38489e00-541d-11eb-8843-a966a400d9ca.png)


An example of a function to open the files is shown below:

<code>
  
    def open_VolcNet_file(file_path, defo_sources):
      """A file to open a single VolcNet file and extrast the deformation source into a one hot encoded numpy array, 
      and the deforamtion location as a n_ifg x 4 array.  
      Ifgs are masked arrays, in m, with up as positive.  

      Inputs:
          file_path | string or Path | path to fie
          defo_sources | list of strings | names of deformation sources, should the same as the names used in VolcNet

      Returns:
          X | r4 masked array | ifgs, as above. ? x y x x n_channels  
          Y_class | r2 array | class labels, ? x n_classes
          Y_loc | r2 array | locations of signals, ? x 4 (as x,y, width, heigh)

      History:
          2020_01_11 | MEG | Written
      """
      import pickle
      import numpy as np
      import numpy.ma as ma

      # 0: Open the file    
      with open(file_path, 'rb') as f:                                                      # open the real data file
          ifgs = pickle.load(f)                                                                # this is a masked array of ifgs
          ifgs_dates = pickle.load(f)                                                          # list of strings, YYYYMMDD that ifgs span
          pixel_lons = pickle.load(f)                                                          # numpy array of lons of lower left pixel of ifgs
          pixel_lats = pickle.load(f)                                                          # numpy array of lats of lower left pixel of ifgs
          all_labels = pickle.load(f)                                                          # list of dicts of labels associated with the data.  e.g. deformation type etc.  
      f.close()        

      # 1: Initiate arrays
      n_ifgs = ifgs.shape[0]                                                                      # get number of ifgs in file
      X = ifgs                                                                                    # soft copy to rename
      Y_class = np.zeros((n_ifgs, len(defo_sources)))                                             # initiate
      Y_loc = np.zeros((n_ifgs, 4))                                                               # initaite

      # 2: Convert the deformation classes to a one hot encoded array and the locations to an array
      for n_ifg in range(n_ifgs):                                                                 # loop through the ifgs
          current_defo_source = all_labels[n_ifg]['deformation_source']                           # get the current ifgs deformation type/class label
          arg_n = defo_sources.index(current_defo_source)                                         # get which number in the defo_sources list this is
          Y_class[n_ifg, arg_n] = 1                                                               # write it into the correct position to make a one hot encoded list.  
          Y_loc[n_ifg, :] = all_labels[n_ifg]['deformation_location']                             # get the location of deformation.  

      return X, Y_class, Y_loc

</code>
