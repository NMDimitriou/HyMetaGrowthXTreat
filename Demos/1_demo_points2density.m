%% Compute density from spatial distributions of cells
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

clear; clc; close all;
addpath('../Points2Density/')

dx=15; % the approximate size of the cell in μm
% 2500μm x 2500μm x 917μm the approximate size of the space
sz = ceil([2500/dx,2500/dx,917/dx]); % scale to the cell size

% Let's create some cells
%center = [2500/2, 2500/2, 917/2];
%coords1 = randi(500,[1000,1]) + center(1);
%coords2 = randi(500,[1000,1]) + center(2);
%coords3 = randi(500,[1000,1]) + round(center(3));
%plot3(coords1,coords2,coords3,'.')
%writematrix([coords1 coords2 coords3],'test_coordinates.txt')


% For the interpolation
x=linspace(0,2.5,sz(1));
z=linspace(0,0.917,sz(3));
[X,Y,Z] = ndgrid(x,x,z);
szq = [480,480,176]; % the number of grid points in the simulation domain
xq=linspace(0,2.5,szq(1));
zq=linspace(0,0.917,szq(3));
[Xq,Yq,Zq] = ndgrid(xq,xq,zq);

%disp(['Initializing density estmation with dsxy=' num2str(dsxy) ', dz=' num2str(dsz)]);
disp(['Initializing density estmation with dxyz=' num2str(dx) 'μm']);
%% Import coordinates
file='test_coordinates.txt';
disp(['Importing coordinates from ' file])
coord = readmatrix(file);
coord = coord./dx; % scale such that cells become points
%shift all coordinates by 1 to remove zeros
coord(:,1:2) = coord(:,1:2)+1;

plot3(coord(:,1),coord(:,2),coord(:,3),'.')
%% Convert centroids to density
disp('Converting points to density...')
[X1,X2,X3] = meshgrid(1:sz(1),1:sz(2),1:sz(3));
d          = 3;
grid       = reshape([X1(:),X2(:),X3(:)],sz(1)*sz(2)*sz(3),d);

dmt = akde(coord,grid);
denscell = reshape(dmt,size(X1));
denscell = denscell.*(dx^3); % convert PDF to Probability
disp('Interpolating to much the simulation grid ...' )
F = griddedInterpolant(X,Y,Z,denscell,'linear');
PV= F(Xq,Yq,Zq);
disp(['min PV = ' num2str(min(PV(:)))])
disp(['min Prob = ' num2str(min(denscell(:)))])
disp(['max PV = ' num2str(max(PV(:)))])
disp(['max Prob = ' num2str(max(denscell(:)))])

%% Plot
disp('Plotting density...')
figure;
view(3)
cmap=flipud(winter(50))';
depth = 0.2;
isovalue3 = depth*(max(PV(:))-min(PV(:)))+min(PV(:));
surf3     = isosurface(xq,xq,zq,PV,isovalue3);
p3        = patch(surf3);
isonormals(xq,xq,zq,PV,p3);
set(p3,'FaceColor',cmap(:,end-3),'EdgeColor','none','FaceAlpha',0.4);   

%% Save matrices to binary files

disp('Saving...')
stat = mkdir('IC'); 
fileid = fopen(['IC/corr_dens.bin'],'w'); 
fwrite(fileid,PV,'double');
fclose(fileid);

disp('Finished!')
%//////