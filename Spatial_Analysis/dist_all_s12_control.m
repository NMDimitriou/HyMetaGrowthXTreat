%% Calculate and plot distances of nuclei for all time points, in each sample
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020
clear; clc; close all;
% the parent directory that contains the subdirectories of the time-points
file1e  = 'res_coord_scaled/';
file2e  = 'res_coord_series_2_scaled/';
file1s  = 'run14_S_LOG_LOCAL_DENSITY_HYBRID_NO_R_NO_DEATH_adhes-5_phenchange0_set12_best_params/';
%file2s  = 'run2_TREAT_HYBRID_NO_R_adhes-5_phenchange0/';
file_s  = [file1s 'scale_coord/'];
%file_s2 = [file2s 'scale_coord/'];


% the list of the time-point subdirectories that contain the data
tp     = {'D0' 'D2'  'D5' 'D7' 'D9' 'D12' 'D14'}; %'D0' 'D2' 
time   = [0 2 5 7 9 12 14]; %0 2  

% the name of the data files for the coordinates
% Series 1
group1 = {'A*C*.txt','A*E*.txt','B*E*.txt','B*N*.txt','B*W*.txt','F*W*.txt'};

% Series 2
group2 = {%'A*C*.txt','A*N*.txt','A*S*.txt','A*W*.txt','A*E*.txt',...
          %'B*C*.txt','B*N*.txt','B*S*.txt','B*W*.txt','B*E*.txt',...
          %'C*C*.txt','C*N*.txt','C*S*.txt','C*W*.txt','C*E*.txt',...
          %'DD*C*.txt','DD*N*.txt','DD*S*.txt','DD*W*.txt','DD*E*.txt',...
                    'E*N*.txt','E*S*.txt','E*W*.txt'           ,...
                    'F*N*.txt','F*W*.txt','F*E*.txt'
          };

groupc1 = {'Control_s1_AC','Control_s1_AE','Control_s1_BE',...
    'Control_s1_BN','Control_s1_BW','Control_s1_FW'};

groupcs1={'CA_coord_Control_ACs1*.txt','CA_coord_Control_AEs1*.txt','CA_coord_Control_BEs1*.txt',...
          'CA_coord_Control_BNs1*.txt','CA_coord_Control_BWs1*.txt','CA_coord_Control_FWs1*.txt'};

groupcs2={'CA_coord_Control_ENs2*.txt','CA_coord_Control_ESs2*.txt','CA_coord_Control_EWs2*.txt',...
          'CA_coord_Control_FNs2*.txt','CA_coord_Control_FWs2*.txt','CA_coord_Control_FEs2*.txt'};

prf_gc2 = 'CA_coord_';
groupc2 = {%'Pac_0p5_AC','Pac_0p5_AN','Pac_0p5_AS','Pac_0p5_AW','Pac_0p5_AE',...
    %'Pac_0p05_BC','Pac_0p05_BN','Pac_0p05_BS','Pac_0p05_BW','Pac_0p05_BE',...
    %'Pac_0p005_CC','Pac_0p005_CN','Pac_0p005_CS','Pac_0p005_CW','Pac_0p005_CE',...
    %'Pac_0p0005_DC','Pac_0p0005_DN','Pac_0p0005_DS','Pac_0p0005_DW','Pac_0p0005_DE',...
    'Control_s2_EN', 'Control_s2_ES', 'Control_s2_EW',...
    'Control_s2_FN', 'Control_s2_FW', 'Control_s2_FE' 
    };

%gname1 = {'Sample_1', 'Sample_2', 'Sample_3', 'Sample_4' , 'Sample_5' , 'Sample_6' };
%gname2 = {'Sample_7', 'Sample_8', 'Sample_9', 'Sample_10', 'Sample_11', 'Sample_12' };

%gname_ca1 = {'CA_Sample_1', 'CA_Sample_2', 'CA_Sample_3', 'CA_Sample_4', 'CA_Sample_5', 'CA_Sample_6' };
%gname_ca2 = {'CA_Sample_7', 'CA_Sample_8', 'CA_Sample_9', 'CA_Sample_10', 'CA_Sample_11', 'CA_Sample_12' };

%gname = [gname1, gname2, gname_ca1, gname_ca2];
datgroup={'Experiments','Simulations'};

nsamp1  = [length(group1), length(groupc1)];
nsamp2  = [12, 12];
ntp    = length(tp);

%vname = {'ACs1','AEs1','BEs1','BNs1','BWs1','FWs1',...
%         'ENs2','ESs2','EWs2','FNs2','FWs2','FEs2'};

run plotopt.m

%% Import coordinates
disp('1. Importing coordinates...')
% Series 1-2
samp1    = cell(nsamp1(1),ntp);
sampc1   = cell(nsamp2(1),ntp);
coord    = struct;
count    = struct;
ca_coord = struct;
[samp1,sampc1,coord,count,ca_coord] = import_coord(nsamp1(2),tp,ntp,file1e,...
    file_s,group1,groupc1,append(prf_gc2,groupc1,'*.txt'),append(prf_gc2,groupc1),...
    samp1,sampc1,coord,count,ca_coord,0);
[samp1,sampc1,coord,count,ca_coord] = import_coord(nsamp1(2),tp,ntp,file2e,...
    file_s,group2,groupc2,append(prf_gc2,groupc2,'*.txt'),append(prf_gc2,groupc2),...
    samp1,sampc1,coord,count,ca_coord,6);%(end-5:end)
%
%[samp1,sampc1,coord,count,ca_coord] = import_coord(nsamp2(2),tp,ntp,file2e,...
%    file_s2,group2(1:20),groupc2(1:20),append(prf_gc2,groupc2(1:20),'*.txt'),append(prf_gc2,groupc2(1:20)),...
%    samp1,sampc1,coord,count,ca_coord,12);
                                          

disp('Finished importing coordinates...')

%% Plot histogram of cells across z
disp('2. Plotting histograms...')
run plotopt.m

%plot_hist_z(coord,ca_coord,tp,'hist_z_s12_exp_ca')
%plot_heat_z(coord,ca_coord,time,917,21,[file_s 'heat_z_s12_exp_ca']);

c=0;
a=fieldnames(ca_coord);
b=fieldnames(coord);
tmp=struct;
expt=struct;
dose={'Control','0p5','0p05','0p005','0p0005'};
ds = {'Control','0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$'};
iters=[1 12; 13 17; 18 22; 23 27; 28 32];
for i=1:1
    for j=iters(i,1):iters(i,2)
        expt.(b{j}) = getfield(coord,b{j});
        disp(a{j})
        tmp.(a{j})  = getfield(ca_coord,a{j});
        disp('----')
    end
    plot_heat_z(expt,tmp,time,917,21,[file_s 'heat_z_s12_exp_ca_' dose{i}], ds{i});
    c = c+5;
    clear tmp expt
end


disp('4. Finished plotting histograms...')

%% Calculate distances and save them
disp('4. Calculating and saving distances...')
gnam = [gname1;gname2];
%parpool(nsamp(1)); %memory is a thing here
for i=1:(nsamp1(1)+nsamp2(2))
    for k=1:ntp
   
        calc_ind_knd(   coord.(gnam{i}).(tp{k}), samp1{i,k}.name)
        calc_ind_knd(ca_coord.(gnam{i}).(tp{k}),sampc1{i,k}.name)
    
    end
end
disp('Finished calculating and saving distances!')
%delete(gcp('nocreate'))

%% Import distance distributions
disp('5. Importing distance distributions...')

IND = struct;
KND = struct;

[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Control_s1_','KND_Control_s1_',group1,groupc1,tp,IND,KND);
[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Control_s2_','KND_Control_s2_',group2,groupc2,tp,IND,KND); %(end-5:end)
%[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Pac_*','KND_Pac_*',group2(1:5),groupc2(1:5),tp,IND,KND);
%[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Pac_*','KND_Pac_*',group2(6:10),groupc2(6:10),tp,IND,KND);
%[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Pac_*','KND_Pac_*',group2(11:15),groupc2(11:15),tp,IND,KND);
%[IND,KND] = import_dist(['Distances/Distributions/'],'IND_Pac_*','KND_Pac_*',group2(16:20),groupc2(16:20),tp,IND,KND);
%
[IND,KND] = import_dist([file1s 'Distances/Distributions/' ],'IND_','KND_',append(prf_gc2,groupc1,'*.txt'),append(prf_gc2,groupc1),tp,IND,KND);
[IND,KND] = import_dist([file1s 'Distances/Distributions/' ],'IND_','KND_',append(prf_gc2,groupc2,'*.txt'),append(prf_gc2,groupc2),tp,IND,KND); %(end-5:end)
%[IND,KND] = import_dist([file2s 'Distances/Distributions/' ],'IND_','KND_',append(prf_gc2,groupc2(1:20),'*.txt'),append(prf_gc2,groupc2(1:20)),tp,IND,KND);

%{
group=groupcs1;
disp('Importing distances...')
for i=1:length(group)
    
    sampINDist{i}=dir(['Distances/Distributions/IND_' group{i} ]);
    sampKNDistCA{i}=dir(['Distances/KNDist/' 'KNDist_' group{i} '.bin']);
    namesINDist  = {sampINDist{i}.name};
    namesINDist  = natsort(namesINDist);
    namesKNDist  = {sampKNDistCA{i}.name};
    namesKNDist  = natsort(namesKNDist);
    
    for j=1:length(tp)
       
     %   fileID=fopen(['Distances/Distributions/' namesINDist{j}],'r');
     %   INDist.(gname{i}).(tp{j})=fread(fileID,'double');
     %   fclose(fileID);
        
        fileID=fopen(['Distances/KNDist/' namesKNDist{j}],'r');
        KNDistCA.(gname{i}).(tp{j})=fread(fileID,'double');
        fclose(fileID);
        
    end
end
disp('Finished')
%}


%% Plot inter-nucleic, and nearest neighbor distance
disp('4. Plotting IN-, NN- distances...')

for i=1:nsamp1(1)+nsamp2
    h=figure('Name',['IND_KND_' gname1{i}],'Position', [200, 200, 1200, 900]);
    alphaval=0.8;
    cmap=winter(7);
    for k=1:ntp
        subplot(1,2,1);
        
        hax=gca;
        hold on
        xi1=IND.(gname1{i}).(tp{k})(2,:);
        f1 =IND.(gname1{i}).(tp{k})(1,:);
        h1=plot(xi1./1000,f1,'LineWidth',4,'Color',cmap(k,:));
        h1.Color(4) = alphaval; 
        xlabel('IN Distances $(mm)$')
        ylabel('PDF') 
             
        subplot(1,2,2);
        hbx=gca;
        hold on
        xi2=KND.(gname1{i}).(tp{k})(2,:);
        f2 =KND.(gname1{i}).(tp{k})(1,:); 
        h2=plot(xi2./1000,f2,'LineWidth',4,'Color',cmap(k,:));
        h2.Color(4) = alphaval; 
        xlabel('NN Distances $(mm)$') 
        colormap(cmap)
        caxis([0 14])
        colorbar(hbx,'Ticks',[0,14],'TickLabels',{'D0' 'D14'},'TickLabelInterpreter','latex'); %
    end
    hold off
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(h,['test_IN_NN_' gname1{i} '.png'],'-dpng','-r350')
end
disp('Finished plotting IN-, NN- distances!')

%% Plot all the samples together
% IND experiments
gname = [groupc2 groupc1];
datgroup = {'Control'}; %'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$', 
nsamp1 = [12]; %5 5 5 5 
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','IND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.65]); %
   
    for j=1:length(gname)
        xi1(:,j)=IND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =IND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    %axis([0 3000 0 inf])
    axis tight
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Experiments','xoff',0.,'yoff',.1);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'IN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.015 .04 .065 .09 .115],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file2s 'IN_all_exp_group12.png'],'-dpng','-r350')

%% IND simulations
gname = [groupc2 groupc1];
gname = append(prf_gc2,gname);
datgroup = {'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$', 'Control'};
nsamp1 = [5 5 5 5 12];
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','IND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.65]); %
   
    for j=1:length(gname)
        xi1(:,j)=IND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =IND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    %axis([0 3000 0 inf])
    axis tight
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Simulations','xoff',0.,'yoff',.08);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'IN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.015 .04 .065 .09 .115],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file2s 'IN_all_sim_group12.png'],'-dpng','-r350')



%% KND experiments
gname = [groupc2 groupc1];
datgroup = {'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$', 'Control'};
nsamp1 = [5 5 5 5 12];
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','KND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.65]); %
   
    for j=1:length(gname)
        xi1(:,j)=KND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =KND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    axis([0 70 0 inf])
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Experiments','xoff',0.,'yoff',.1);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'NN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.015 .04 .065 .09 .115],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file2s 'NN_all_exp_group12.png'],'-dpng','-r350')

%% KND simulations
gname = [groupc2 groupc1];
gname = append(prf_gc2,gname);
datgroup = {'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$', 'Control'};
nsamp1 = [5 5 5 5 12];
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','KND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.65]); %
   
    for j=1:length(gname)
        xi1(:,j)=KND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =KND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    axis([0 70 0 inf])
    %axis tight
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Simulations','xoff',0.,'yoff',.08);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'NN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.015 .04 .065 .09 .115],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file2s 'NN_all_sim_group12.png'],'-dpng','-r350')

%%
% KND experiments
gname = groupc2;
datgroup = {'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$'};
nsamp1 = [5 5 5 5];
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','KND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.5]);
   
    for j=1:length(gname)
        xi1(:,j)=KND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =KND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    %axis([0 3000 0 inf])
    %axis tight
    axis([0 149 0 inf])
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Experiments','xoff',0.,'yoff',.1);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'NN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.02 .06 .1 .14],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file 'KN_all_exp_group2.png'],'-dpng','-r350')

%% KND simulations
gname = append(prf_gc2,groupc2);
datgroup = {'0.5 $\mu M$','0.05 $\mu M$','0.005 $\mu M$','0.0005 $\mu M$'};
nsamp1 = [5 5 5 5];
clear f1
run plotopt.m
overlap = 0.6;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','KND','Position', [100, 100, 1200, 500]);
for i=1:length(tp)
    
    subplot(1,5,i,'Position',[0.1+0.9*(i-1)/6 0.2 0.12 0.5]);
   
    for j=1:length(gname)
        xi1(:,j)=KND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =KND.(gname{j}).(tp{i})(1,:);
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1));% + nsamp1(1:ii-1));
        else
            m=0;
        end
        
        for j=1:(nsamp1(ii))%+nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1) 
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    title(tp{i})
    %axis([0 3000 0 inf])
    axis([0 149 0 inf])
    %axis tight
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005'});%,'Control'
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit('Simulations','xoff',0.,'yoff',.1);
%set(hh.th,'yoff',.01);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',20);
ylh=ylabel(hh.ah,'NN Distances');
set(ylh, 'Visible', 'On','FontSize',30)
colormap(cmap)
c=colorbar('Ticks',[.02 .06 .1 .14],'TickLabels',... %0.1290/9 0.1290/4
    datgroup,...
    'TickLabelInterpreter','latex',...
    'Position',[0.85 0.15 0.02 0.5]); %0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file 'KN_all_sim_group2.png'],'-dpng','-r350')

%% IND
run plotopt.m
gname = [groupc1 groupc2 append(prf_gc2,groupc1) append(prf_gc2,groupc2)]; % 
overlap = 0.4;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','IND','Position', [100, 100, 1000, 900]);
for i=1:length(tp)
    
    subplot(2,4,i);
   
    for j=1:length(gname)
        xi1(:,j)=IND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =IND.(gname{j}).(tp{i})(1,:);
        f1(:,j) = f1(:,j)/max(f1(:,j));
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1) + nsamp1(1:ii-1));
        else
            m=0;
        end
        for j=1:(nsamp1(ii) + nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1)
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    disp('===========')
    title(tp{i})
    axis([0 3000 0 inf])
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005','Control'});
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',33);
ylh=ylabel(hh.ah,'Inter-Nucleic Distances');
set(ylh, 'Visible', 'On','FontSize',33)
colormap(cmap)
c=colorbar('Ticks',[0.1290/9 0.1290/4],'TickLabels',...
    {'Experiments','Simulations'},...
    'TickLabelInterpreter','latex',...
    'Position',[0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566]);
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file_s 'IN_all_exp_sim_group12.png'],'-dpng','-r350')

%% KND
run plotopt.m
gname = [groupc1 groupc2 append(prf_gc2,groupc1) append(prf_gc2,groupc2)]; % 
overlap = 0.4;
mn = linspace(0,1,length(gname));
cmap=winter(length(datgroup));

h=figure('Name','KND','Position', [100, 100, 1000, 900]);
for i=1:length(tp)
    
    subplot(2,4,i);
   
    for j=1:length(gname)
        xi1(:,j)=KND.(gname{j}).(tp{i})(2,:);
        f1(:,j) =KND.(gname{j}).(tp{i})(1,:);
        f1(:,j) = f1(:,j)/max(f1(:,j));
    end
    f = cumsum(max(f1,[],1))*(1-overlap);
    
    hold on
    for ii=length(datgroup):-1:1
        if(ii>1)
            m=sum(nsamp1(1:ii-1) + nsamp1(1:ii-1));
        else
            m=0;
        end
        for j=1:(nsamp1(ii) + nsamp1(ii))
            patch([xi1(:,j+m)',-min(xi1(:,j+m))],[f1(:,j+m)+f(ii);f(ii)],mn(ii),...
                 'EdgeColor','none','FaceAlpha',0.05,'FaceColor',cmap(ii,:));

            plot(xi1(:,j+m),f1(:,j+m)+f(ii),'Color',cmap(ii,:),'LineWidth',1)
            disp(['j+m = ' num2str(j+m)])
        end
        disp('---------')
    end
    disp('===========')
    title(tp{i})
    axis([0 149 0 inf])
    yticks(f)
    yticklabels({' '})
    %if(i==1 || i==5)
    %    yticklabels({'Pac 0.5','Pac 0.05','Pac 0.005','Pac 0.0005','Control'});
    %else
    %    yticklabels({' '})
    %end
    hold off
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
set(xlh, 'Visible', 'On','FontSize',33);
ylh=ylabel(hh.ah,'Nearest-Neighbour Distances');
set(ylh, 'Visible', 'On','FontSize',33)
colormap(cmap)
c=colorbar('Ticks',[0.1290/9 0.1290/4],'TickLabels',...
    {'Experiments','Simulations'},...
    'TickLabelInterpreter','latex',...
    'Position',[0.5422+1.3*0.1566  0.11  0.02  0.11+1.25*0.1566]);
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%      'Location','southeastoutside','NumColumns',2);
%newPosition = [0.4 0.2 0.2 0.1];
%newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
%title(hh.ah,'Inter-Nuclei distance')
print(h,[file_s 'KN_all_exp_sim_group12.png'],'-dpng','-r350')

%% Estimate variations between distributions
% import distances
disp('5. Importing distances...')
for i=1:(nsamp1(1)+nsamp2(2))
    for k=1:ntp
    
    sampINDist{i,k}=dir(['Distances/INDist/' 'INDist_' samp1{i,k}.name '.bin' ]);
    sampKNDist{i,k}=dir(['Distances/KNDist/' 'KNDist_' samp1{i,k}.name '.bin']);
         
    fileID=fopen([sampINDist{i,k}.folder '/' sampINDist{i,k}.name],'r');
    INDist.(gname1{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);
        
    fileID=fopen([sampKNDist{i,k}.folder '/' sampKNDist{i,k}.name],'r');
    KNDist.(gname1{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);
      
    end
end
for i=1:(nsamp1(1)+nsamp2(2))
    for k=1:ntp
    
    sampINDistCA{i,k}=dir(['Distances/INDist/' 'INDist_' sampc1{i,k}.name '.bin' ]);
    sampKNDistCA{i,k}=dir(['Distances/KNDist/' 'KNDist_' sampc1{i,k}.name '.bin']);
         
    fileID=fopen([sampINDistCA{i,k}.folder '/' sampINDistCA{i,k}.name],'r');
    INDist.(gname1{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);
        
    fileID=fopen([sampKNDistCA{i,k}.folder '/' sampKNDistCA{i,k}.name],'r');
    KNDist.(gname_ca1{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);
      
    end
end

disp('Finished importing distances!')
clear sampINDist sampKNDist sampINDistCA sampKNDistCA

%% calculate the cosine similarity between distributions
disp('6. Estimating similarity...')
%parpool(length(group));
for i=1:nsamp
   
    calc_dist_var(INDist.(gname1{i}),KNDist.(gname1{i}),gname1{i});
    
end
%delete(gcp('nocreate'))
disp('Finished estimating similarity!')

%% Import variations between distributions 
disp('7. Importing variations between distributions...')


for i=1:nsamp1(1)%+nsamp2(2)
 
    varIND.(groupc1{i})=readmatrix([file1s 'Distances/Variability/Var_INDist_Sim_v_Exp' groupc1{i} '.txt']);
    varIND.(groupc2{i})=readmatrix([file1s 'Distances/Variability/Var_INDist_Sim_v_Exp' groupc2{i} '.txt']);

    varKND.(groupc1{i})=readmatrix([file1s 'Distances/Variability/Var_KNDist_Sim_v_Exp' groupc1{i} '.txt']);
    varKND.(groupc2{i})=readmatrix([file1s 'Distances/Variability/Var_KNDist_Sim_v_Exp' groupc2{i} '.txt']);

end

%% Sim vs exp

% find total mean and std
vind = cell2mat(struct2cell(varIND));
vknd = cell2mat(struct2cell(varKND));

meanIND_sim = mean(vind(:,4))
meanKND_sim = mean(vknd(:,4))

stdIND_sim = std(vind(:,4))
stdKND_sim = std(vknd(:,4))


%%
varINDt=struct2array(varIND);
varINDt=vertcat(varINDt(:,1:3),varINDt(:,4:6),varINDt(:,7:9),varINDt(:,10:12),varINDt(:,13:15),varINDt(:,16:18));

varKNDt=struct2array(varKND);
varKNDt=vertcat(varKNDt(:,1:3),varKNDt(:,4:6),varKNDt(:,7:9),varKNDt(:,10:12),varKNDt(:,13:15),varKNDt(:,16:18));

varINDmean  = mean(varINDt(:,3)); 
varKNDmean  = mean(varKNDt(:,3));

varINDstd  = std(varINDt(:,3));
varKNDstd  = std(varKNDt(:,3));

disp('Finished importing variations between distributions!')
%% Plot similarity for sim vs exp

mat_ind = reshape(vind(:,4),[7 nsamp1(1)+nsamp1(2)]);
mat_knd = reshape(vknd(:,4),[7 nsamp1(1)+nsamp1(2)]);

mean_mat_ind = mean(mat_ind,2);
mean_mat_knd = mean(mat_knd,2);

std_mat_ind = std(mat_ind,0,2);
std_mat_knd = std(mat_knd,0,2);

run plotopt.m

h=figure('Name','Variations','Position', [100, 100, 1000, 800]);   
hold on
errorbar(time,mean_mat_ind,std_mat_ind,'LineWidth',5,'DisplayName','Inter-Nucleic Distances')
errorbar(time,mean_mat_knd,std_mat_knd,'LineWidth',5,'DisplayName','Nearest-Neighbour Distances')
axis([0 15 0 1.1])
xlabel('time (days)','FontSize',44)
ylabel('cosine similarity','FontSize',44)
leg=legend('Location','southwest','FontSize',40);
set(leg,'color','none');
hold off

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,[file1s 'var_IN_KN_sim_vs_exp.png'],'-dpng','-r350')

%{
%% Plot similarity for different time-points
disp('8. Plotting similarity...')

vIND_dt = [];
vKND_dt = [];
dt  = [];

for i=1:nsamp1(1)+nsamp2(2)
    
    %i-tp difference
    subin{i}=find(abs(varINDt(:,1)-varINDt(:,2))==i);
    subkn{i}=find(abs(varKNDt(:,1)-varKNDt(:,2))==i);
    len = length(subkn{i});

    vIND_dt=[vIND_dt;    varINDt(subin{i},3)]; 
    vKND_dt=[vKND_dt;    varKNDt(subkn{i},3)];
    dt     =[dt; i*ones(len,1)]; 
    
end

h=figure('Name','Variations','Position', [100, 100, 1200, 900]);   
hold on
violinplot(vIND_dt,dt,'ShowMean',true,'ViolinAlpha',0.5);
xlabel('$\Delta t_{i,j}$')
ylabel('cosine similarity')
%axis([-inf inf 0.9 1])
hold off
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'test_violin_IN_s1.png','-dpng','-r350')

%
h=figure('Name','Variations','Position', [100, 100, 1200, 900]);   
hold on
violinplot(vKND_dt,dt,'ShowMean',true,'ViolinAlpha',0.5);
xlabel('$\Delta t_{i,j}$')
ylabel('cosine similarity')
%axis([-inf inf 0.9 1])
hold off
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'test_violin_NN_s1.png','-dpng','-r350')

disp('Finished plotting similarity!')

%% Plot all the samples together
disp('9. Plotting all distance samples...')
% IND
h=figure('Name','IND','Position', [100, 100, 1200, 900]);
for i=1:ntp
    
    subplot(3,3,i);
    %subaxis(3,3,i, 'Spacing', 0.08, 'Padding', 0.015, 'Margin', 0.08);
    hax=gca;
    hold on
    for j=1:nsamp
        xi1=IND.(gname1{j}).(tp{i})(2,:);
        f1 =IND.(gname1{j}).(tp{i})(1,:);
        h1=plot(xi1./1000,f1,'LineWidth',4);
        box off
    end
    title(tp{i})
    axis([-inf inf 0 10e-04])
    hax.FontSize=24;
    hold off
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(mm)$','FontSize',30);
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'PDF','FontSize',30);
set(ylh, 'Visible', 'On');
hL=legend(gname1,'Interpreter','latex','FontSize',24,...
      'Location','southeastoutside','NumColumns',2);
newPosition = [0.4 0.15 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
print(h,['test_IND_group1.png'],'-dpng','-r350')
%title(hh.ah,'Inter-Nucleic distance')

% KND
h=figure('Name','KND','Position', [100, 100, 1200, 900]);
for i=1:ntp
    
    subplot(3,3,i);
    %subaxis(3,3,i, 'Spacing', 0.08, 'Padding', 0.015, 'Margin', 0.08);
    hax=gca;
    hold on
    for j=1:nsamp
        xi1=KND.(gname1{j}).(tp{i})(2,:);
        f1 =KND.(gname1{j}).(tp{i})(1,:);
        h1=plot(xi1./1000,f1,'LineWidth',4);
        hax.FontSize=24;
    end
    title(tp{i})
    axis([-inf 200/1000 0 0.1])
    box off
    hold off
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(mm)$','FontSize',30);
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'PDF','FontSize',30);
set(ylh, 'Visible', 'On');
hL=legend(gname1,'Interpreter','latex','FontSize',24,...
      'Location','southeastoutside','NumColumns',2);
newPosition = [0.4 0.15 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
print(h,['test_KND_group1.png'],'-dpng','-r350')
disp('Finished plotting all distance samples...')
%--------------------------------------------------------------------------

%}