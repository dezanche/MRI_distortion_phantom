# MRI_distortion_phantom

This modular phantom uses solid signal sources to measure distortion in magnetic resonance imaging (MRI) systems. Its large volume and light weight make it practical for applications such as radiation therapy where accurate patient contours are critical for treatment planning. Please note that image distortion can depend on variables such as sequence parameters and magnetic field (B<sub>0</sub>) homogeneity.
This repository contains detailed information in addition to the scientific abstract submitted to the [2022 ISMRM conference](https://www.ismrm.org/22m/).
  
## Contents
[CAD](./tree/main/CAD) contains CAD files describing the phantom as built.\
[Images](./tree/main/Images) contains screenshots and photographs.\
[Support_files](./tree/main/Support_files) contains software and other files used for processing data obtained with the phantom.\
[Construction](#construction) describes how the phantom is designed and built.\
Images acquired using the phantom can be processed to analyze the distortion using both [3D Slicer](#data-processing-using-3d-slicer) and [Matlab](#data-processing-using-matlab).

## Construction
The phantom is made using off-the-shelf parts as much as possible. The [CAD](./tree/main/CAD) drawing was made in [FreeCAD](https://www.freecad.org/) and includes all elements of the phantom as built. Parameters in the table included in the file can be changed to adapt the model to different dimensions.\
The proton signal is generated by silicone rubber beads (9 mm ⌀ spheres with a 2 mm ⌀ axial hole for mounting), attached to a support structure using M2 nylon screws and nuts.
The support structure consists of stacked polypropylene sheets, 3.3 mm thick and with 4 mm ⌀ holes at regular intervals (25.4 mm “peg board”). The sheets are joined using nylon threaded rods and separated by stand-offs 22.1 mm in height, resulting in a 3-dimensional grid with isotropic pitch.
The present configuration uses 11 sheets with 11⨉18 holes, for a possible 2178 sample locations.
Slotted side panels can be added to increase flexural, torsional, and shear stiffness.
Three rows from the edges of the boards were populated with 1386 beads to allow a 25 mm ⌀ container filled with ~12 mL of petroleum jelly to be placed at the centre of the phantom for automatic flip angle and frequency calibration. More positions can be populated if desired.\
![scanner](https://github.com/dezanche/MRI_distortion_phantom/blob/main/Images/Photos/20211028_103933.jpg)
**ADD Rendering and close-up photo OF PHANTOM**

## Data Processing Using 3D Slicer 

After scanning the volume occupied by the phantom, export the volume in DICOM format and analyze in [3D Slicer](https://www.slicer.org/) as follows.

### 1. Data Import
- DICOM browser module: select series and load
- Volume rendering module: display preset > MRI MIP > adjust "shift" to see cloud of dots
![MIP](https://github.com/dezanche/MRI_distortion_phantom/blob/main/Images/Screenshots/2021-10-29-Scene_801_MIP.png)

### 2. Segment
- Segment editor module: add segment > threshold > set threshold* > apply > show 3D \
(adjust low range of threshold while viewing RGY screens: start high and lower it to see beads but minimize noise and artefacts) \
then clean up using
- Segment editor module: islands > Split islands to segments> set min size > apply \
(e.g., min size 50 voxels)
- (optional) clean up manually by excluding segments corresponding to, e.g., noise, artefacts, central sample, etc.
![segmentation](https://github.com/dezanche/MRI_distortion_phantom/blob/main/Images/Screenshots/2021-10-29-Scene_801_segment.png)

### 3. Calculate Centroid of Each Segment
- following the methods in (https://discourse.slicer.org/t/centroid-determination/3541) copy and paste the following code into Slicer's Python Interactor:
```python
segmentationNode = getNode("Segmentation")  # or whatever name was given to the segmentation
# Compute centroids
import SegmentStatistics
segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
segStatLogic.getParameterNode().SetParameter("Segmentation", segmentationNode.GetID())
segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.centroid_ras.enabled", str(True))
segStatLogic.computeStatistics()
stats = segStatLogic.getStatistics()
# Place a markup point in each centroid
markupsNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
markupsNode.CreateDefaultDisplayNodes()
for segmentId in stats["SegmentIDs"]:
  centroid_ras = stats[segmentId,"LabelmapSegmentStatisticsPlugin.centroid_ras"]
  segmentName = segmentationNode.GetSegmentation().GetSegment(segmentId).GetName()
  markupsNode.AddFiducialFromArray(centroid_ras, segmentName)
```
  
- lock markups corresponding to centroids so they can't be accidentally moved with the mouse:\
Markups module > Control Points > Interaction in Views > lock
- uncheck *Point Labels* under *Advanced* to view the cloud of fiducials without labels cluttering the scene

### 4. Import Known Grid Positions
- create an FCSV file containing the XYZ coordinates of all beads in the phantom using `bead_coordinates.m' or other means (see the [Markups module documentation](https://www.slicer.org/wiki/Documentation/4.10/Modules/Markups) and https://discourse.slicer.org/t/convert-csv-files-to-be-loadable-as-fiducials-or-as-sequence-of-fiducials-json/18320/2)
- load the FCSV file directly from UI using “Add data” dialog
- lock markups corresponding to beads so they can't be accidentally moved with the mouse:\
Markups module > Control Points > Interaction in Views > lock

### 5. Manual Alignment
If some rotation and translation are needed to align the 2 sets of fiducials proceed as follows:
- Data module > right-click “Transform” column (on the beads coordinates) > Create new transform
- right-click Transform > Edit properties...
- adjust translation and rotation sliders to get the 2 sets aligned
- to make the changes permanent: Data module > Transform hierarchy tab > right click on the transform > "harden transform"
- save transformed coordinates to FCSV file if desired
![fiducials](https://github.com/dezanche/MRI_distortion_phantom/blob/main/Images/Screenshots/2021-11-01-Scene_801_fiducials.png)

### 6. Calculate Displacements
The following code (paste into Slicer's Python Interactor) calculates displacements between each measured centroid and the corresponding known coordinates, discards pairs more than `maxDistance` apart, and writes the results in a table.
```python
import scipy
import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
print('finding nearest neighbour pairs...\n')

maxDistance = 20;   #distance beyond which points are ignored [mm]

knownPointsNode = getNode('bead_coordinates');   # known points 3D Slicer node
queryPointsNode = getNode('MarkupsFiducial');   # query points 3D Slicer node

knownPoints = slicer.util.arrayFromMarkupsControlPoints(knownPointsNode);     # fiducials of known grid points
queryPoints = slicer.util.arrayFromMarkupsControlPoints(queryPointsNode);    # fiducials of measured centroids

# create the k-dimensional tree from the known points:
tree = KDTree(knownPoints);

# query the tree with the query points using upper bound maxDistance:
distances, indices = tree.query(queryPoints, k=1, distance_upper_bound=maxDistance);

print('minimum distance =', min(distances), ' \n')

# select pairs that are not assigned inf distance
nearIndices = indices[np.isfinite(distances)];
nearQueryPoints = queryPoints[np.isfinite(distances)];
reorderedKnownPoints = np.take(knownPoints,nearIndices,axis = 0);
distanceMatrix = cdist(nearQueryPoints,reorderedKnownPoints);
reorderedDistances = distanceMatrix.diagonal();

radialCoordinate = np.linalg.norm(reorderedKnownPoints, axis=1);

print('maximum distance =', max(reorderedDistances), ' \n')

# write result into table; see https://www.slicer.org/wiki/Documentation/4.10/ScriptRepository#Create_histogram_plot_of_a_volume
tableNode=slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
tableNode.SetName('Coordinates and Displacements')
updateTableFromArray(tableNode, np.hstack((nearQueryPoints, reorderedKnownPoints, np.vstack((radialCoordinate, reorderedDistances)).T, np.absolute(reorderedKnownPoints-nearQueryPoints))),["actual X","actual Y","actual Z","grid X","grid Y","grid Z","rho","displacement","|dX|","|dY|","|dZ|"]);

# launch pop-up table
tableView=slicer.qMRMLTableView()   
tableView.setMRMLTableNode(tableNode)
tableView.show()
```

### 7. Plot Results
The data in the table can be saved to a TSV file and analyzed using other software, or plotted directly in 3D Slicer as follows.
- Plots module > Charts tab > Create New Plotchart, Create New Plotseries
- Series tab > select Data Series, Plot Type > Scatter, Input Table > Coordinates and Displacements
- choose X and Y Axis columns (e.g., the radial coordinate *rho* and *displacement* magnitude, respectively, in the example below)
- Line Style > none
- if not visible, display the plot window by clicking the eye in the Charts tab\
or main menu: View > Layout > 3D Table
![scatter](https://github.com/dezanche/MRI_distortion_phantom/blob/main/Images/Screenshots/2021-11-01-Scene_801_scatter.png)

## Data Processing Using Matlab
The acquired data can also be processed in Matlab using [Distortion_Analysis.m](./blob/main/Support_files/Distortion_Analysis.m).
Please note that this code will not run in [GNU Octave](https://www.gnu.org/software/octave/) because a couple of critical functions used in the script are not available in Octave toolboxes.
