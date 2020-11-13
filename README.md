# trajectory_visualization: 


This is a [Shiny App](https://www.shinyapps.io) to visualize trajectory data of Intrinsically Disordered Proteins/Regions (IDPs/IDRs). You have the option to visualize an IDR attached to a folded domain (FD) or by itself. 



## Usage

Clone or download the repository to whatever directory you please. Open the file `shiny_app.R` in RStudio and click `Run App` in the top right corner. A window will open prompting you to enter the following parametres: 

* `IDR Start Residue` - Residue where IDR begins
* `IDR End Residue` - Residue where IDR ends
* `FD Start Residue` - Residue where FD starts
* `FD End Residue` - Residue where FD ends
* `Num Clusters` - Number of clusters to split each frame in (see Methods)
* `RMSD Resolution` - Number of equally spaced alpha-carbon atoms to calculate the RMSD (see Methods)
* `Frame Stride` - Spacing between frames. For instance, if this parameter equals 10, we will use every 10th frame. The application can comfortably handle around 1000 frames to visualize, and up to 3000 if you really want to push it. Adjust this parameter depending on the size of your simulation
* `Display Mean Trajectory` - Display each trajectory from each frame or display the mean trajectory (see Methods)
* `File to PDB` - Location to PDB file 
* `File to Traj` - Location to trajectory file (expects either xtc or dcd format)
* `Folder to Traj` - Location to a parent folder underneath which are other folders containing trajectory data. For instance, if you have data organized in the following format, passing in the path to the folder `traj_data`, will recursively search the directory and get all files of the format `xtc` and `dcd`. 

###### Data Organization Tree
    traj_data
    |
    +----> 1
           |
           +----> traj.xtc
    |
    +----> 2
           |
           +----> traj.xtc
    |
    +----> 3
           |
           +----> traj.xtc

In order for the visualization to render, every parameter needs to be entered with a valid value and *either* `File to Traj` or `Folder to Traj`



## Packages 
For R, the following packages are needed:

* `shiny` 
* `shinyFiles` 
* `dplyr`
* `plotly`
* `reshape2`

For Python, the following packages are needed:

* `numpy` 
* `scipy` 
* `pandas`
* `sklearn`
* `mdtraj`
* `Bio.PDB`


## Methods

There's two main components to generating the input to visualize the trajectory data: 

### Trajectory Processing

If dealing with multiple trajectories, all trajectories are first [joined] [2] into a single trajectory. Then, each frame in the trajectory is [aligned] [3] to the first frame. If the simulation consists of an IDR and FD, only the coordinates of the FD are used to align each frame. If a FD is present, the aligned coordinates of the FD from the first frame are used to visualize it. 

### Trajectory Clustering

To aid in the visualization of trajectories with a large number of frames, each frame is treated as independent instance and [hierarchically clustered] [4] according to its [Root-mean-square deviation] [5] (RMSD). Specifically, the RMSD is only calculated with respect to alpha-carbon atoms. The user has the option to specify the number of equally spaced alpha-carbon atoms they want to use to calculate the RMSD using the argument `rmsd_ca_resolution`. For instance, if an IDP is 30 residues long and `rmsd_ca_resolution = 10`, the alpha-carbon of the 3rd, 6th, 9th, etc. residues will be used to calculate the RMSD. This option is provided because the coordinates of residues in a polymer are correlated to a certain extent (as quantitivaley defined by the [persistence length] [6]), so one may get 'better' clusters as they tune that parameter. [Ando et. al] [7] cite a persistence length of a C-terminal tail around 1nm (roughly 3 residues). Granted, this is for a specific system, but provides a good starting point for playing around with this parameter. 

### Trajectory Visualization

In addition to visualizing each individual frame, the user also has the option to visualize the mean trajectory in each cluster. This is accomplished by averaging the `<x,y,z>` coordinates of each residue across all frames for a given cluster. Also plotted alongside the mean trajectory is the mean trajectory ± 2 SD. As of now, a mean trajectory is only plotted for clusters with more than 50 instances.  


## Examples

Below are images of the GFP protein with an IDR attached at the N-terminal region. Residues on the FD are color coded based on their charge (red = negative, grey = neutral, blue = positive). Opening the image in a new tab will display higher resolution versions of them. 

[Figure 1] [8] -- All frames displayed for all cluster

[Figure 2] [9] -- All frames displayed for an individual cluster

[Figure 3] [10] -- Mean trajectory displayed for all clusters







[2]: https://mdtraj.org/1.9.4/api/generated/mdtraj.join.html?highlight=join#mdtraj.join
[3]: https://mdtraj.org/1.9.4/api/generated/mdtraj.Trajectory.superpose.html?highlight=superpose#mdtraj.Trajectory.superpose 
[4]: https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering
[5]: https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
[6]: https://en.wikipedia.org/wiki/Persistence_length
[7]: https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/cphc.200800210
[8]: https://imgur.com/jf1e1hp
[9]: https://imgur.com/xLYrXV9
[10]: https://imgur.com/pStTl1X



