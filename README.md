# Dislocation3DAnalysis
4D Analysis of DFXM Data using Matlab

To run the analysis in this repository, the general workflow is to run the following scripts in sequence of each other:
(1) "Segmentation.m"
(2) "plotSegmentedDislocationStructures.m"
(3) "BoundaryAnalysis.m"

In each script, you will need to follow the instructions in the comments, making sure to clearly define the folder where the data resides (always called "dataPath") and the folder in which you wish to run save all the saved output files (always called "writePath"). 

STEP 1: IMAGE SEGMENTATION
In the first script, you will need to input a path to the full array of all your data files from the 2D z,phi or z,chi scans, and the begin by inputting your experimental parameters from the ESRF ID06-HXM (e.g. photon energy, magnification, etc). You will then need to manually uncomment the different beginning code blocks to resolve the micrometers per pixel from your experiment (to calibrate the position axes) and the images that define the weak-beam condition for each z-layer. 

Based on the memory usage in your computer and your scan size, you may need to scale the image arrays. There are scaling factors when you load your files into Matlab, and another set before you begin the computationally-intensive fiber texture methods that are required to connect the disconnected pixels that define the dislocations into linear segments for further analysis.

Once the images are saved into memory and the fiber texturing is complete, you will also need to manually define your intensity threshold to define what sets a "dislocation pixel" vs a "background pixel". This can be done by changing the first element of the 2-element array called "clim", which is initially defined in the code as "0.1". Finally, the Segmentation script will save all the relevant variables onto your computer for use in the next script.

STEP 2: PLOTTING MOVIE AND IMAGES OF 3D DISLOCATION STRUCTURES
The script to plot the dislocation structures is relatively unsupervised. Once the output from your Segmentation method saves, it should run without issue once you update the writePath, dataPath, and filename to match the file saved by Segmentation. Following this, you will need to update the quiver text to match the crystallographic orientations you used in your experiment, so that your crystallographic assignments match the legend in your plots. Finally, you may want to change the number of frames in the video, the frame rate, or the quality it saves with - all of which are in the final block of code in that script. The Matlab figure should save along with the 3D spinning movie - though the Boolean variables at the beginning allow you to turn on/off the save features for these functions.

STEP 3: BOUNDARY ASSIGNMENT VIA DIMENSIONAL REDUCTION
To assign the boundaries, you must begin the "BoundaryAnalysis" script by defining the x,y,z bounds that define each boundary you are interested in, and the colors that you would like each of them plotted in (in order, RGB vectors for the colors). The script will then open a parallel loop (set the number of nodes if that matters to you) and will define the point as presnt or not present in neach boundary array you have defined. Note that you will need to manually set the number of if/elseif statements in the second block of code to ensure that the parfor loop can run without carrying too much memory to each node. 

You will then need to define your coordinate transformation matrices (detailed in the Supplementary Information) to account for the crystal orientation vectors and the diffraction angles in the third block of code. The script will then use the MSAC algorithm to find the boundary plane normals and plot the results of the probabilistic method. You will need to run this a number of times to evaluate what the convergence of the method is. Note that you will also need to set the value of "maxDistance" based on the noise of your measurement. That variable will define the distance out of the plane beyond which the method will declare any point to be an outlier and exclude it from the fit. A value for this that is too low will result in poor boundary identification, with a mean error that is essentially the same as the maxDistance value; as you increase that value, you should eventually find that the MSAC mean error stops changing based on the maxDistance size. That value is the one appropriate for an accurate result.

In the fourth block of code, you will then need to parse the "Boundary" structure variable that is defined above and identify the values for the hkl and the error that are most commonly found when you run the method a reasonable number of times (~30x). Note that you will need to round the numbers and decide what maximum value of hkl is reasonable to approximate as infinity (i.e. zero), and this is best done with guidance from the plots compiled in this block of code.

Finally, in the fifth block of code, the pointClouds are dimensionally reduced into their 2D planes and 1D line vectors. This block of code is largely manual, so you will need to ensure that you go through it carefully to ensure that the solutions are accurate for your data. I left a few commented lines of code that plot additional traces that can be helpful to understand the math and transformations in the code. Once you identify the correct boundary angles, you will again need to input them manually, as here is no automated angle identification method.

Once you finish running "BoundaryAnalysis", you will need to parse through the structure variable "Boundary", as all of the results of the algorithms will be saved as fields in that variable.