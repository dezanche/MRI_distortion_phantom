# MRI_distortion_phantom

These files describe a large-volume, lightweight, modular phantom with solid signal sources used to measure distortion in magnetic resonance images.

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
- following the methods in (https://discourse.slicer.org/t/centroid-determination/3541) paste the following code into Slicer's Python interactor:
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
- create FCSV file containing the XYZ coordinates of all beads in the phantom using `bead_coordinates.m' or other means (see the [Markups module documentation](https://www.slicer.org/wiki/Documentation/4.10/Modules/Markups) and https://discourse.slicer.org/t/convert-csv-files-to-be-loadable-as-fiducials-or-as-sequence-of-fiducials-json/18320/2)
- load the FCSV file directly from UI using “Add data” dialog
- lock markups corresponding to beads so they can't be accidentally moved with the mouse:\
Markups > Control Points > Interaction in Views > lock
