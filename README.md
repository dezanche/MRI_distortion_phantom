# MRI_distortion_phantom

These files describe a large-volume, lightweight, modular phantom with solid signal sources used to measure distortion in magnetic resonance imaging (MRI) systems. Please note that image distortion can depend on variables such as sequence parameters and magnetic field (B_0) homogeneity.

## Data Processing Using 3D Slicer

After scanning the volume occupied by the phantom, export the volume in DICOM format and analyze in [3D Slicer](https://www.slicer.org/) as follows.

### 1. Data Import
- DICOM browser module: select series and load
- Volume rendering module: display preset > MRI MIP > adjust "shift" to see cloud of dots

### 2. Segment
- Segment editor module: add segment > threshold > set threshold* > apply > show 3D \
(adjust low range of threshold while viewing RGY screens: start high and lower it to see beads but minimize noise and artefacts) \
then clean up using
- Segment editor module: islands > Split islands to segments> set min size > apply \
(min size 50 voxels)
- (optional) clean up manually by excluding segments corresponding to, e.g., noise, artefacts, central sample, etc.

### 3. Calculate Centroid of Each Segment
- following the methods in (https://discourse.slicer.org/t/centroid-determination/3541) paste the following code into Slicer's Python Interactor:
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
- choose X and Y Axis columns (e.g., *rho* and *displacement*, respectively)
- Line Style > none
- if not visible, display the plot window by clicking the eye in the Charts tab or main menu: View > Layout > 3D Table
