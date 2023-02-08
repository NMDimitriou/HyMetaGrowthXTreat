%% 3D plot of centroids
clear; clc;



tp    = {'D5' 'D7' 'D9' 'D12' 'D14'}; %'D0' 'D2' 
ts    = {'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14'}; %'D0' 'D1' 'D2' 'D3' 'D4' 
time  = [5 7 9 12 14];% 0 2 

% Series 1
%group1 = {'AD*C.tif.txt','AD*E.tif.txt','BD*E.tif.txt','BD*N.tif.txt','BD*W.tif.txt','FD*W.tif.txt'};
prf_s  = 'num_cells_';
prf_c  = 'ca_res_';
%prf_ctr1 = 'Control_s1_';
%groups1= {'AC','AE','BE','BN','BW','FW'}; %num_cells_Control_s1_FW.txt

% CA series 1,2
%
% Series 2
%prf_ctr2 = 'Control_s2_';
prf_ctr2  = 'Pac_';
group2 = {'A*C*.txt','A*N*.txt','A*S*.txt','A*W*.txt','A*E*.txt',...
          'B*C*.txt','B*N*.txt','B*S*.txt','B*W*.txt','B*E*.txt',...
          'C*C*.txt','C*N*.txt','C*S*.txt','C*W*.txt','C*E*.txt',...
          'D*C*.txt','D*N*.txt','D*S*.txt','D*W*.txt','D*E*.txt',...
          %          'E*N*.txt','E*S*.txt','E*W*.txt'           ,...
          %          'F*E*.txt','F*N*.txt','F*W*.txt'
         };

groups2 = {'0p5_AC','0p5_AN','0p5_AS','0p5_AW','0p5_AE',...
    '0p05_BC','0p05_BN','0p05_BS','0p05_BW','0p05_BE',...
    '0p005_CC','0p005_CN','0p005_CS','0p005_CW','0p005_CE',...
    '0p0005_DC','0p0005_DN','0p0005_DS','0p0005_DW','0p0005_DE'
    %'EN','ES','EW','FE','FN','FW'
    };                
          

run plotopt.m
group_sim = 'run3_RADIAL_TREAT_HYBRID_NO_R_adhes-5_phenchange0/'; %'run2_TREAT_HYBRID_NO_R_adhes-5_phenchange0/';
file1e  = 'res_coord_scaled/';
file2e  = 'res_coord_series_2_scaled/';
%file_s  = 'res_coord_sim_series_12/run3_HYBRID_NO_R_adhes-5_phenchange0/';
iters=1;
ntp     = length(tp);

%%

for i=1:length(group1)
   
    for j=1:length(tp)
       
        samp1{i,j}=dir(['res_coord_scaled/' tp{j} '/' group1{i}]);
        coord1.([prf_ctr1 groups1{i}]).(tp{j})=readmatrix([samp1{i,j}.folder '/' samp1{i,j}.name]);
        count1.([prf_ctr1 groups1{i}]).(tp{j})=length(coord1.([prf_ctr1 groups1{i}]).(tp{j})(:,1));

    end
    for k=1:iters
        samps1{i}=dir([group_sim prf_s prf_ctr1 groups1{i} '.txt']); %'_iter_' num2str(k)
        sampc1{i}=dir([group_sim prf_c prf_ctr1 groups1{i} '.csv']);% '_iter_' num2str(k)
        ca_count1.([prf_ctr1 groups1{i}]) = readmatrix([samps1{i}.folder '/'  samps1{i}.name]);%{k}
        ca_coord1.([prf_ctr1 groups1{i}]) = readmatrix([sampc1{i}.folder '/'  sampc1{i}.name]);
    end
    
end
%%
for i=1:length(group2)  
    for j=1:length(tp)
       
        samp2{i,j}=dir(['res_coord_series_2_scaled/' tp{j} '/' group2{i}]);
        coord2.([prf_ctr2 groups2{i}]).(tp{j})=readmatrix([samp2{i,j}.folder '/' samp2{i,j}.name]);
        count2.([prf_ctr2 groups2{i}]).(tp{j})=length(coord2.([prf_ctr2 groups2{i}]).(tp{j})(:,1));

    end
    for k=1:iters
        samps2{i}=dir([group_sim prf_s prf_ctr2 groups2{i}  '.txt']); %'_iter_' num2str(k) 
        sampc2{i}=dir([group_sim prf_c prf_ctr2 groups2{i}  '.csv']); %'_iter_' num2str(k) 
        ca_count2.([prf_ctr2 groups2{i}]) = readmatrix([samps2{i}.folder '/'  samps2{i}.name]);
        ca_coord2.([prf_ctr2 groups2{i}]) = readmatrix([sampc2{i}.folder '/'  sampc2{i}.name]);%{k}
    end
end

%%
% Split time-points of CA data
%ca_coord1_tsp = split_tp(append(prf_ctr1,groups1), append(prf_ctr1,groups1), ca_coord1, ts, iters);
ca_coord2_tsp = split_tp(append(prf_ctr2,groups2), append(prf_ctr2,groups2), ca_coord2, ts, iters);

%%

% Reformat coordinates
sz=[480,480,176];
h=2500/480; %scale factor

for i=1:length(ts)
    mkdir([group_sim 'scale_coord/' ts{i} '/']);
end

%

%ca_coord1_tsp = save_sim_split_coord(group_sim,append(prf_ctr1,groups1),ca_coord1_tsp,append(prf_ctr1,groups1),ts,h,sz,iters);
ca_coord2_tsp = save_sim_split_coord(group_sim,append(prf_ctr2,groups2),ca_coord2_tsp,append(prf_ctr2,groups2),ts,h,sz,iters);


%% Merge controls and count all
gc1=struct2array(count1);
gc2=struct2array(count2);
gc1=struct2table(gc1);
gc2=struct2table(gc2);
gc1=table2array(gc1);
gc2=table2array(gc2);

%gc_control = gc1;
%gc_control = gc2;
gc_control = [gc1; gc2];
%gc_control = [gc1; gc2(21:26, :)];
%gc_pac_0p5 = gc2(1:5,:);
%gc_pac_0p05 = gc2(6:10,:);
%gc_pac_0p005 = gc2(11:15,:);
%gc_pac_0p0005 = gc2(16:20,:);

gcm_control   = mean(gc_control);
gcstd_control = std(gc_control);

%gcm_pac_0p5 = mean(gc_pac_0p5);
%gcm_pac_0p05 = mean(gc_pac_0p05);
%gcm_pac_0p005 = mean(gc_pac_0p005);
%gcm_pac_0p0005 = mean(gc_pac_0p0005);

%gcstd_pac_0p5 = std(gc_pac_0p5);
%gcstd_pac_0p05 = std(gc_pac_0p05);
%gcstd_pac_0p005 = std(gc_pac_0p005);
%gcstd_pac_0p0005 = std(gc_pac_0p0005);

% find std of ca simulated data

ts = (time(1):0.1:time(end));

% interpolate accoring to ts
%%
ca_count_table = [];
ca_count_table = merge_counts(ca_count_table,append(prf_ctr1,groups1),...
    ca_count1,append(prf_ctr1,groups1),ts,0);

%
ca_count_table = merge_counts(ca_count_table,append(prf_ctr2,groups2),...
    ca_count2,append(prf_ctr2,groups2),ts,6); %was 6

ca_count_mean_ctr = mean(ca_count_table,2);
ca_count_std_ctr  = std(ca_count_table,0,2);

%ca_count_0p5    = ca_count_table(:,1:5);
%ca_count_0p05   = ca_count_table(:,6:10);
%ca_count_0p005  = ca_count_table(:,11:15);
%ca_count_0p0005 = ca_count_table(:,16:20);
%

%ca_count_mean_0p5= mean(ca_count_0p5,2);
%ca_count_mean_0p05= mean(ca_count_0p05,2);
%ca_count_mean_0p005= mean(ca_count_0p005,2);
%ca_count_mean_0p0005= mean(ca_count_0p0005,2);

%ca_count_std_0p5 = std(ca_count_0p5,0,2);
%ca_count_std_0p05 = std(ca_count_0p05,0,2);
%ca_count_std_0p005 = std(ca_count_0p005,0,2);
%ca_count_std_0p0005 = std(ca_count_0p0005,0,2);


%% Plot growth for controls
f1=figure('Name','Growth','Position', [100, 100, 600, 500]);
hold on
h1=plot(ts,ca_count_mean_ctr,'Color','#00308F','LineWidth',5,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean_ctr'-ca_count_std_ctr';ca_count_mean_ctr'+ca_count_std_ctr'],'#00308F');
h2=errorbar(time,gcm_control,gcstd_control,'.-','Color','#D3212D','LineWidth',5,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','northwest','FontSize',22)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time (days)','Interpreter','latex','FontSize',34)
ylabel('Nuclei count','Interpreter','latex','FontSize',34)
title('Control')
%legend('Location','NorthWest')
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[group_sim 'nuclei_count_controls_s12_no_death.png'],'-r350','-dpng')


%% Plot growth for 0.5\mu M
f1=figure('Name','Growth','Position', [100, 100, 600, 500]);
hold on
h1=plot(ts,ca_count_mean_0p5,'Color','#00308F','LineWidth',5,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean_0p5'-ca_count_std_0p5';ca_count_mean_0p5'+ca_count_std_0p5'],'#00308F');
h2=errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','Color','#D3212D','LineWidth',5,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','southwest','FontSize',22)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time (day)','Interpreter','latex','FontSize',34)
ylabel('Nuclei count','Interpreter','latex','FontSize',34)
title('0.5 $\mu m$')
%legend('Location','NorthWest')
axis([5 14 0 inf])
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[group_sim 'nuclei_count_s2_0p5.png'],'-r350','-dpng')

%% 0.05
f1=figure('Name','Growth','Position', [100, 100, 600, 500]);
hold on
h1=plot(ts,ca_count_mean_0p05,'Color','#00308F','LineWidth',5,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean_0p05'-ca_count_std_0p05';ca_count_mean_0p05'+ca_count_std_0p05'],'#00308F');
h2=errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','Color','#D3212D','LineWidth',5,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','southwest','FontSize',22)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time (day)','Interpreter','latex','FontSize',34)
ylabel('Nuclei count','Interpreter','latex','FontSize',34)
title('0.05 $\mu m$')
%legend('Location','NorthWest')
axis([5 14 0 inf])
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[group_sim 'nuclei_count_s2_0p05.png'],'-r350','-dpng')
%% 0.005
f1=figure('Name','Growth','Position', [100, 100, 600, 500]);
hold on
h1=plot(ts,ca_count_mean_0p005,'Color','#00308F','LineWidth',5,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean_0p005'-ca_count_std_0p005';ca_count_mean_0p005'+ca_count_std_0p005'],'#00308F');
h2=errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','Color','#D3212D','LineWidth',5,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','southwest','FontSize',22)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time (day)','Interpreter','latex','FontSize',34)
ylabel('Nuclei count','Interpreter','latex','FontSize',34)
title('0.005 $\mu m$')
%legend('Location','NorthWest')
axis([5 14 0 inf])
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[group_sim 'nuclei_count_s2_0p005.png'],'-r350','-dpng')


%% 0.0005
f1=figure('Name','Growth','Position', [100, 100, 600, 500]);
hold on
h1=plot(ts,ca_count_mean_0p0005,'Color','#00308F','LineWidth',5,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean_0p0005'-ca_count_std_0p0005';ca_count_mean_0p0005'+ca_count_std_0p0005'],'#00308F');
h2=errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','Color','#D3212D','LineWidth',5,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','northwest','FontSize',22)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time (day)','Interpreter','latex','FontSize',34)
ylabel('Nuclei count','Interpreter','latex','FontSize',34)
title('0.0005 $\mu m$')
%legend('Location','NorthWest')
axis([5 14 0 inf])
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[group_sim 'nuclei_count_s2_0p0005.png'],'-r350','-dpng')
%% Plot centroids Series 1-2 experiment
plot_centr(group1,gname1,coord1,tp,'centroids');
plot_centr(group2,gname2,coord2,tp,'centroids');


clf; close all

%% Plot centroids Series 1-2 CA
plot_centr_CA(groupc1,gname1,ca_coord1_tsp,tp,'centroids_CAs1');
plot_centr_CA(groupc2,gname2,ca_coord2_tsp,tp,'centroids_CAs2');


clf; close all

%% Plot CA and experiment together
run plotopt.m
plot_centr_CA_exp(append(prf_ctr1,groups1),ca_coord1_tsp,coord1,append(prf_ctr1,groups1),tp,[group_sim 'centroids_CA_exp_s1'])
plot_centr_CA_exp(append(prf_ctr2,groups2),ca_coord2_tsp,coord2,append(prf_ctr2,groups2),tp,[group_sim 'centroids_CA_exp_s2'])


clf; close all