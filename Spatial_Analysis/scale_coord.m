%% Scale coordinates

clear; clc;
%% Series 1
%{
tp    = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
group = {'A*C*.txt','A*E*.txt','B*E*.txt','B*N*.txt','B*W*.txt','F*W*.txt'};
gname = {'AC','AE','BE','BN','BW','FW'};

% import
for i=1:length(group)
   
    for j=1:length(tp)
       
        samp{i,j}=dir(['res_coord_scaled/' tp{j} '/' group{i}]);
        coord.(gname{i}).(tp{j})=readmatrix([samp{i,j}.folder '/' samp{i,j}.name]);

    end
end

% scale
scalexy = 0.4023; %p ixels/microns
scalez  = 0.05;   % slices/microns

for i=1:length(group)
   
    for j=1:length(tp)
       
        coord.(gname{i}).(tp{j})(:,1:2)=coord.(gname{i}).(tp{j})(:,1:2)/scalexy;
        coord.(gname{i}).(tp{j})(:,3  )=coord.(gname{i}).(tp{j})(:,3  )/scalez ;
        writematrix(coord.(gname{i}).(tp{j}),[samp{i,j}.folder '/' samp{i,j}.name]);

    end
end
%}

%% Series 2
tp    = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
group = {'A*C*.txt','A*N*.txt','A*S*.txt','A*W*.txt','A*E*.txt',...
         'B*C*.txt','B*N*.txt','B*S*.txt','B*W*.txt','B*E*.txt',...
         'C*C*.txt','C*N*.txt','C*S*.txt','C*W*.txt','C*E*.txt',...
         'D*C*.txt','D*N*.txt','D*S*.txt','D*W*.txt','D*E*.txt',...
                    'E*N*.txt','E*S*.txt','E*W*.txt'           ,...
                    'F*N*.txt',           'F*W*.txt','F*E*.txt'};
gname = {'AC','AN','AS','AW','AE',...
         'BC','BN','BS','BW','BE',...
         'CC','CN','CS','CW','CE',...
         'DC','DN','DS','DW','DE',...
              'EN','ES','EW'     ,...
              'FN',     'FW','FE'};

% import
for i=1:length(group)
   
    for j=1:length(tp)
       
        samp{i,j}=dir(['res_coord_series_2/' tp{j} '/' group{i}]);
        coord.(gname{i}).(tp{j})=readmatrix([samp{i,j}.folder '/' samp{i,j}.name]);

    end
end

% scale
scalexy = 0.4023; % pixels/micron
scalez  = 0.0667;   % slices/micron

for i=1:length(group)
   
    for j=1:length(tp)
       
        wrt{i,j}=dir(['res_coord_series_2_scaled/' tp{j} '/' group{i}]);
        coord.(gname{i}).(tp{j})(:,1:2)=coord.(gname{i}).(tp{j})(:,1:2)/scalexy;
        coord.(gname{i}).(tp{j})(:,3  )=coord.(gname{i}).(tp{j})(:,3  )/scalez ;
        writematrix(coord.(gname{i}).(tp{j}),[wrt{i,j}.folder '/' samp{i,j}.name]);

    end
end