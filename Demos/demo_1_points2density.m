%% Compute density from spatial distributions of cells
% Author: Nikolaos M. Dimitriou,
% McGill University, 2020

clear; clc; close all;
addpath('../Points2Density/')

dx=15; % the approximate size of the cell in μm
% 2500μm x 2500μm x 917μm the approximate size of the space
sz = ceil([2500/dx,2500/dx,917/dx]); % scale to the cell size

% Let's create some cells
%
center = [2500/2, 2500/2, 917/2];
coords1 = randi([-400,400],[1000,1]) + center(1);
coords2 = randi([-400,400],[1000,1]) + center(2);
coords3 = randi([-300,300],[1000,1]) + round(center(3));
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D0.txt')

coords1 = randi([-400,400],[900,1]) + center(1);
coords2 = randi([-400,400],[900,1]) + center(2);
coords3 = randi([-300,300],[900,1]) + round(center(3)/1.1);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D2.txt')

coords1 = randi([-600,600],[900,1]) + center(1);
coords2 = randi([-600,600],[900,1]) + center(2);
coords3 = randi([-150,200],[900,1]) + round(center(3)/3);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D5.txt')

coords1 = randi([-600,600],[800,1]) + center(1);
coords2 = randi([-600,600],[800,1]) + center(2);
coords3 = randi([-150,200],[800,1]) + round(center(3)/3);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D7.txt')

coords1 = randi([-700,700],[500,1]) + center(1);
coords2 = randi([-700,700],[500,1]) + center(2);
coords3 = randi([-114,200],[500,1]) + round(center(3)/4);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D9.txt')

coords1 = randi([-700,700],[500,1]) + center(1);
coords2 = randi([-700,700],[500,1]) + center(2);
coords3 = randi([-90,180],[500,1]) + round(center(3)/5);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D12.txt')

coords1 = randi([-700,700],[300,1]) + center(1);
coords2 = randi([-700,700],[300,1]) + center(2);
coords3 = randi([-90,120],[300,1]) + round(center(3)/5);
figure;
plot3(coords1,coords2,coords3,'.')
writematrix([coords1 coords2 coords3],'test_coordinates_D14.txt')


% For the interpolation
x=linspace(0,2.5,sz(1));
z=linspace(0,0.917,sz(3));
[X,Y,Z] = ndgrid(x,x,z);
szq = [480,480,176]; % the number of grid points in the simulation domain
xq=linspace(0,2.5,szq(1));
zq=linspace(0,0.917,szq(3));
[Xq,Yq,Zq] = ndgrid(xq,xq,zq);
dsxy = 2500/480; % downscale parameter for xy-plane
dsz  = 917/176;  % downscale parameter for z-dimension

%disp(['Initializing density estmation with dsxy=' num2str(dsxy) ', dz=' num2str(dsz)]);
disp(['Initializing density estmation with dxyz=' num2str(dx) 'μm']);
%% Import coordinates
%file='test_coordinates.txt';
fnames = 'test_coordinates_';
day_names = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
ndays = length(day_names);
coord={};
for i=1:ndays
    disp(['Importing coordinates from ' day_names{i}])
    coord.(day_names{i}) = readmatrix([fnames day_names{i} '.txt']);
    coord.(day_names{i})(:,1:2) = ceil(coord.(day_names{i})(:,1:2) ./dsxy); % scale such that cells become points
    coord.(day_names{i})(:,3  ) = ceil(coord.(day_names{i})(:,3  ) ./dsz);
    %shift all coordinates by 1 to remove zeros
    coord.(day_names{i})(:,1:2) = coord.(day_names{i})(:,1:2)+1;

    figure;
    plot3(coord.(day_names{i})(:,1),coord.(day_names{i})(:,2),coord.(day_names{i})(:,3),'.')
end
%% Convert centroids to density
disp('Converting points to density...')
[X1,X2,X3] = meshgrid(1:szq(1),1:szq(2),1:szq(3));
d          = 3;
grid       = reshape([X1(:),X2(:),X3(:)],szq(1)*szq(2)*szq(3),d);
PV = {};
for i=1:ndays
    dmt = akde(coord.(day_names{i}),grid);
    PV.(day_names{i}) = dmt*length(coord.(day_names{i}))*(15^3)/(5.21^3);
    disp(['min PV = ' num2str(min(PV.(day_names{i})(:)))])
    disp(['max PV = ' num2str(max(PV.(day_names{i})(:)))])
    PV.(day_names{i}) = reshape(PV.(day_names{i}), [480 480 176]);
end


%% Plot
disp('Plotting density...')

cmap=flipud(winter(50))';
depth = 0.2;

for i=1:ndays
    disp([' -> ' day_names{i}])
    figure;
    view(3)
    isovalue3 = depth*(max(PV.(day_names{i})(:))-min(PV.(day_names{i})(:)))+min(PV.(day_names{i})(:));
    surf3     = isosurface(xq,xq,zq,PV.(day_names{i}),isovalue3);
    p3        = patch(surf3);
    isonormals(xq,xq,zq,PV.(day_names{i}),p3);
    set(p3,'FaceColor',cmap(:,end-3),'EdgeColor','none','FaceAlpha',0.4);
    xlim([0 2.5])
    ylim([0 2.5])
    zlim([0 0.917])
end

%% Save matrices to binary files

disp('Saving...')
stat = mkdir('IC');
for i=1:ndays
    fileid = fopen(['IC/new_dens_' day_names{i} '.raw'],'w');
    fwrite(fileid,PV.(day_names{i}),'double');
    fclose(fileid);
end
disp('Finished!')
%//////
