%%
clear; clc;


fd1 = '../../IC/series1/';
fd2 = '../../IC/series2/';
tp  = ['D2','D5','D7','D9','D12','D14'];
tim = [2 5 7 9 12 14];
prf1= 'Control_s1_';
prf2= 'Control_s2_';
prf = 'new_dens_';


dat1 = ['AC','AE','BE','BN','BW','FW'];
dat2 = ['EN','ES','EW','FN','FW','FE'];

%%

% Control 1
for i=1:length(dat1)
    disp(['Importing control ' dat1(i)]);
    for j=1:length(tim)
   
        tmp = dir([fd1 prf1 dat1(i) '/' prf prf1 dat1(i) '_' tp(j) '.raw']);
        fileid = fopen([tmp.folder '/' tmp.name]);
        A = fread(fileid,'double');
        fclose(fileid);
    
        dat{i,j} = A;
    
        clear A;
    end
end

% Control 2
for i=1:length(dat2)
    disp(['Importing control ' dat2(i)]);
    for j=1:length(tim)
   
        tmp = dir([fd2 prf2 dat2(i) '/' prf prf2 dat2(i) '_' tp(j) '.raw']);
        fileid = fopen([tmp.folder '/' tmp.name]);
        A = fread(fileid,'double');
        fclose(fileid);
    
        dat{i+length(dat1),j} = A;
    
        clear A;
    end
end

dat

% Simulations
for i=1:length(dat1)
   
    disp(['Importing sim ' dat1(i)]);
    for j=1:length(tim)
        
        tmp = dir([prf1 dat1(i) '/a_double_' tp(j) '*.raw']);
        fileid = fopen([tmp.folder '/' tmp.name]);
        A = fread(fileid,'double');
        fclose(fileid);
    
        sim{i,j} = A;
    
        clear A;
        
    end
end


for i=1:length(dat2)
   
    disp(['Importing sim ' dat2(i)]);
    for j=1:length(tim)
        
        tmp = dir([prf2 dat2(i) '/a_double_' tp(j) '*.raw']);
        fileid = fopen([tmp.folder '/' tmp.name]);
        A = fread(fileid,'double');
        fclose(fileid);
    
        sim{i+length(dat1),j} = A;
    
        clear A;
        
    end
end

sim


for i=1:length(dat1)+length(dat2)
    
    for j=1:length(tim)

        cdc(i,j) = continuous_dice_coefficient(dat{i,j},sim{i,j});

    end
end

save cdc.mat cdc
