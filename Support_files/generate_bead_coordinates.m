% generates fiducials file (.fcsv) to import bead coordinates into 3D slicer
% see https://www.slicer.org/wiki/Modules:Fiducials-Documentation-3.6
% © 2021 Nicola De Zanche

% coordinate system is rotated to be compatible with 3D Slicer:
% X --> R
% Z --> A
% Y --> S

output_file_name = 'bead_coordinates';

% same variable names used in FreeCAD drawing
xPoints = 11; % # of points along x
yPoints = 16; % less than actual available holes
zPoints = 11;

xPitch = 25.4; % distance [mm] between points along x
yPitch = 25.4;
zPitch = 25.4;

% rectangular exclusion limits to define a box 3 rows deep
xLimit = 60; %[mm]
yLimit = 120;
zLimit = 130;  % top and bottom of box are not filled in

% generate complete 3D grid
xCoords = xPitch*([0:xPoints-1]-(xPoints-1)/2);
yCoords = yPitch*([0:yPoints-1]-(yPoints-1)/2);
zCoords = zPitch*([0:zPoints-1]-(zPoints-1)/2);
[Xall,Yall,Zall] = meshgrid(xCoords,yCoords,zCoords);

% select points in a rectangular shell
indices = find((abs(Xall(:))>xLimit)|(abs(Yall(:))>yLimit)|(abs(Zall(:))>zLimit));
XYZ=[Xall(indices),Yall(indices),Zall(indices)];
total_beads = size(indices,1)

% plot for verification
figure(1)
plot3(Xall(:),Yall(:),Zall(:),'k.');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

figure(2)
plot3(XYZ(:,1),XYZ(:,3),XYZ(:,2),'k.');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

% text to prepend
text = {'# name = bead coordinates'
['# numPoints = ' num2str(total_beads)]
'# symbolScale = 5.5'
'# symbolType = 11'
'# visibility = 1'
'# textScale = 12.5'
'# color = 0.4,1,1'
'# selectedColor = 0.807843,0.560784,1'
'# opacity = 1'
'# ambient = 0'
'# diffuse = 1'
'# specular = 0'
'# power = 1'
'# locked = 1'
'# columns = label,x,y,z,sel,vis'}

fid = fopen([output_file_name '.fcsv'],'w');
fprintf(fid,'%s\n',text{:});  % write file header
fprintf(fid,'\n');
% write coordinates to file
dlmwrite(fid,[[1:total_beads]',XYZ(:,1),XYZ(:,3),XYZ(:,2),zeros(size(indices)),ones(size(indices))],'append','on','delimiter',',');

fclose(fid);