% analyze distortion measured using silicone bead grid phantom
% Copyright ©2021 Mariko Gardiner, Nicola De Zanche 

%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <https://www.gnu.org/licenses/>.

%clear;
%warning off;

bead_radius = 4.5;  % radius of silicone beads [mm]
hole_radius = 0.9;  % radius of hole in silicone beads [mm]
read_data = 1;  % for testing: skip reading data if already in workspace

%% initial checks

% Matlab or Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% check for image toolbox
if isOctave
  hasIPT = license('test', 'image');
  error('this script cannot run in Octave because of missing functions!');
else
  hasIPT = license('test', 'image_toolbox');
end

if ~hasIPT
    % User does not have the toolbox installed.
    message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
    reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
    if strcmpi(reply, 'No')
        % User said No, so exit.
        return;
    end
end
clear hasIPT

%% settings
%Define variables that might be user selected. xr yr and zr control
%rotation of the grid around the x y and z axes respectively, and lim sets the
%threshold for converting the grayscale data to black and white. These all
%may need to vary for different MRI scanners and settings.
xr=0;
yr=0;
zr=0;
limit=16;

%%get data from Dicom to point locations
if read_data
  Volume_path = uigetdir;
  if exist(Volume_path,'file')~=7
    message = fprintf('File not found. \n');
    return
  end
  %user selects method of fitting the data, either by rectangular or
  %spherical distance
  s=listdlg('SelectionMode','Single','ListString',{'Rectangular','Spherical'});
  
  %Read in Dicom volume and metadata
  disp(['loading ' Volume_path ' ...']);
  if isOctave
    [c,sp]=(dicomreadVolume(Volume_path));
  else
    [V,sp]=(dicomreadVolume(Volume_path));
    %Reduce to 3D from 4D
    c=squeeze(V);
  end
end
%If the data comes in a different rotation than the sample data, it is 
%necessary to change the rotation here for the data to keep the xyz axes
%consistent
%c=imrotate3(c,90,[0 1 0]);

%% Create a template sphere using the pixel spacings to determine the size.
% N.B. sp.SliceIncrement doesn't exist in the native Matlab function
if isOctave
  % template about 2 diameters wide
  Nx=round(2*bead_radius/sp.PixelSpacing(1));
  Ny=round(2*bead_radius/sp.PixelSpacing(2));
  Nz=round(2*bead_radius/sp.SliceIncrement);
  
  [X,Y,Z] = meshgrid(sp.PixelSpacing(1)*[-Nx:Nx],sp.PixelSpacing(2)*[-Ny:Ny],sp.SliceIncrement*[-Nz:Nz]);
  ball = X.^2+Y.^2+Z.^2 <= bead_radius^2;
  ball = ball.*(X.^2+Y.^2 >= hole_radius^2);    % TO DO: check orientation of hole in XYZ coordinate system
  
else
  
  ball=false(21,21,11);   % inefficient --> fix; assumes slice increment = PixelSpacing?
  for i=1:21
    for j=1:21
      for k=1:11
        if norm([(11-i)*sp.PixelSpacings(1),(11-j)*sp.PixelSpacings(2),(6-k)])<4.5
          if norm([(11-i)*sp.PixelSpacings(1),(11-j)*sp.PixelSpacings(2)])>2
            ball(i,j,k)=1;
          end
        end
      end
    end
  end
end %if

if isOctave
%  convolved = convn(c,ball);
else
  %make the image binary and filter out noise. The function filters by how
  %many connected pixels there are, so it may be necessary to lower the
  %connectivity for lower resolution images with fewer pixels.
  c=bwareaopen(c>limit,100,26);
  %Label each ball and find the centroid
  V=bwconncomp(c,26);
  p=regionprops3(V,'Centroid');
  f=regionprops3(V,'Volume');
  
  % Find the Dice coefficient for each bead, ignoring beads whose volumes
  % would go outside the volume of the image
  DC=zeros(size(f.Volume,1),1);
  for i=1:size(f.Volume)
    bi=(p.Centroid(i,:));
    if bi(2)<(size(c,1)-10)&&bi(2)>11&&bi(1)>11&&bi(1)<(size(c,2)-10)...
      &&bi(3)>6&&bi(3)<(size(c,3)-6)
      DC(i)=dice(ball,c((bi(2)-10):(bi(2)+10),(bi(1)-10):(bi(1)+10),(bi(3)-5):(bi(3)+5)));
    end
  end
end %if

%adjust the points according to the metadata so that the units are
%milimeters and the grid is approximately centered. If the data is stored
%differently than the sample data, it may be necessary to adjust the equations.
p=p.Centroid;
if sp.PatientPositions(1,1)~=sp.PatientPositions(2,1)
xp=p(:,1 )*(sp.PatientPositions(1,1)-sp.PatientPositions(2,1))-sp.PatientPositions(1,1);
else
  xp=p(:,1 )*sp.PixelSpacings(1)+sp.PatientPositions(1,1); 
end
if sp.PatientPositions(1,3)~=sp.PatientPositions(2,3)
yp=p(:,2 )*(sp.PatientPositions(1,2)-sp.PatientPositions(2,2))-sp.PatientPositions(1,2);
else
  yp=p(:,2 )*sp.PixelSpacings(1)-sp.PatientPositions(1,3); 
end
if sp.PatientPositions(1,2)~=sp.PatientPositions(2,2)
zp=p(:,3 )*(sp.PatientPositions(3,2)-sp.PatientPositions(2,2))-sp.PatientPositions(1,2)*0;
else
  zp=p(:,3 )*sp.PixelSpacings(1)-sp.PatientPositions(1,3)*0; 
end

%% Plot the calculated points in a 3D scatterplot in a separate window
figure('Name','3D plot');
scatter3(xp,yp,zp,'k.');
axis equal
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
hold on
%make the information available in two formats for accessibility
p(:,1)=xp;
p(:,2)=yp;
p(:,3)=zp;


%Create a grid to model the expected values. The dimensions of the grid
%are based off measurements of the physical phantom.
Y=-114.3:25.4:114.3;
X=[-127.0:25.4:-76.2;76.2:25.4:127.0];
Z=10:25.3:263 ;
[x0,y0,z0] = meshgrid(X,Y,Z);
Y1=[-190.5:25.4:-139.7;139.7:25.4:190.5];
X1=-127.0:25.4:127.0;
[x1,y1,z1]=meshgrid(X1,Y1,Z);
p1=[x1(:),y1(:),z1(:)];
p0=[x0(:),y0(:),z0(:)];
p2=[p0;p1;0,0,126.5];
%Rotate the grid if necessary. This is done using the variables defined at
%the beginning
p2=p2*[1 0 0; 0 cos(xr) -sin(xr) ; 0 sin(xr) cos(xr)]*[cos(yr) 0 sin(yr)...
    ; 0 1 0 ; -sin(yr) 0  cos(yr)] *[cos(zr) -sin(zr) 0 ; sin(zr) cos(zr) 0 ...
    ; 0 0 1];

%Find the points on the grid closest to each measured point, and create a
%matrix l with these points in order
[k] = dsearchn(p2,p);
i=1;
l=p;
while i<=size(p,1)
    l(i,:,:)=p2(k(i),:,:);
    i=i+1;
end
%Try to center the grid and recalculate the points closest to each points.
%Limit is the distance from the center within which points are considered.
%It must be greater than 10.
limit=100;
dif=[11, 11, 11];
n=0;
disp('centering grid...');
%This loop calibrates the grid until it is sufficiently matched
while abs(dif(1))>5
    %Based off the user selection earlier, the grid is either calibrated by
    %spherical distance from the center or rectangular distance.
    %Rectangular distance is limited only by y because it is the only
    %variable large enough, but for smaller limits it may be necessary to
    %limit the other directions also. Data outside the limit or inside a
    %smaller limit are not considered when centering the data.
switch s
    case 1
        for i=size(p,1):-1:1
            if abs(p2(i,2))<limit && abs(p2(i,2))>10
                dif=dif+[l(i,1)-xp(i),l(i,2)-yp(i),l(i,3)-zp(i)];
                n=n+1;
            end
        end
    case 2
        for i=size(p,1):-1:1
            if sqrt(p2(i,1).^2+p2(i,2).^2+p2(i,3).^2)<limit && abs(p2(i,2))>10
                dif=dif+[l(i,1)-xp(i),l(i,2)-yp(i),l(i,3)-zp(i)];
                n=n+1;
            end
        end
end
dif=dif/n;
X=p2(:,1,:)-dif(1);
Y=p2(:,2,:)-dif(2);
Z=p2(:,3,:)-dif(3);
p2=[X,Y,Z];
i=1;
[k,d] = dsearchn(p2,p);
while i<=size(p,1)
    l(i,:,:)=p2(k(i),:,:);
    i=i+1;
end
end
%Display an arrow starting from each measured point with arrow size
%determined by the distance to the expected value
quiver3(xp,yp,zp,l(:,1,:)-xp,l(:,2,:)-yp,l(:,3,:)-zp,'off','.b')
hold off

%% Create a new window to display the 2D scatterplots
h=figure('Name','Scatterplots');
tiledlayout(h,'flow')
m=mean(d);
fprintf('%s has an average error of %d . \n',Volume_path,m)
%display scatterplots of error and location based on user selection
e=[l(:,1)-xp,l(:,2,:)-yp,l(:,3,:)-zp];
a=sqrt(e(:,1).^2+e(:,2).^2+e(:,3).^2);
b=sqrt(l(:,1).^2+l(:,2).^2+(l(:,3)-126.5).^2);
s=listdlg('ListString',{'xx','xy','xz','yx'...
    ,'yy','yz','zx','zy','zz',':x',':y',':z','x:','y:','z:','::','Dice coefficient'});
if ismember(1,s)
    scatter(nexttile,l(:,1),e(:,1),'.');
    xlabel('X Location (mm)');
    ylabel('X Error (mm)');
end
if ismember(2,s)
    scatter(nexttile,l(:,1),e(:,2),'.');
    xlabel('X Location (mm)');
    ylabel('Y Error (mm)');
end
if ismember(3,s)
    scatter(nexttile,l(:,1),e(:,3),'.');
    xlabel('X Location (mm)');
    ylabel('Z Error (mm)');
end
if ismember(4,s)
    scatter(nexttile,l(:,2),e(:,1),'.');
    xlabel('Y Location(mm)');
    ylabel('X Error (mm)');
end
if ismember(5,s)
    scatter(nexttile,l(:,2),e(:,2),'.');
    xlabel('Y Location (mm)');
    ylabel('Y Error (mm)');
end
if ismember(6,s)
    scatter(nexttile,l(:,2),e(:,3),'.');
    xlabel('Y Location (mm)');
    ylabel('Z Error (mm)');
end
if ismember(7,s)
    scatter(nexttile,l(:,3),e(:,1),'.');
    xlabel('Z Location(mm)');
    ylabel('X Error (mm)');
end
if ismember(8,s)
    scatter(nexttile,l(:,3),e(:,2),'.');
    xlabel('Z Location (mm)');
    ylabel('Y Error (mm)');
end
if ismember(9,s)
    scatter(nexttile,l(:,3),e(:,3),'.');
    xlabel('Z Location(mm)');
    ylabel('Z Error (mm)');
end
if ismember(10,s)
    scatter(nexttile,b,e(:,1),'.');
    xlabel('Distance from Center(mm)');
    ylabel('X Error (mm)');
end
if ismember(11,s)
    scatter(nexttile,b,e(:,2),'.');
    xlabel('Distance from Center(mm)');
    ylabel('Y Error (mm)');
end
if ismember(12,s)
    scatter(nexttile,b,e(:,3),'.');
    xlabel('Distance from Center(mm)');
    ylabel('Z Error (mm)');
end
if ismember(13,s)
    scatter(nexttile,l(:,1),a,'.');
    xlabel('X Location (mm)');
    ylabel('Total Error (mm)');
end
if ismember(14,s)
    scatter(nexttile,l(:,2),a,'.');
    xlabel('Y Location (mm)');
    ylabel('Total Error (mm)');
end
if ismember(15,s)
    scatter(nexttile,l(:,3),a,'.');
    xlabel('Z Location (mm)');
    ylabel('Total Error (mm)');
end
if ismember(16,s)
    scatter(nexttile,b,a,'.');
    xlabel('Distance from Center (mm)');
    ylabel('Total Error (mm)');
end
if ismember(17,s)
scatter(nexttile,b,DC,'.')
xlabel('Distance from Center');
ylabel('Dice-Sorensen Coefficient');
end
%Clear unnecessary variables. Variables can be added or removed freely.
clear x0 y0 z0 x1 y1 z1 p0 p1 Y1 X1 xr yr zr tf n l1 l2 i j k s;