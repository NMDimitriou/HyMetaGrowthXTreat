%% Plots all the summary functions for three-dimensional point patterns.
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

clear; clc; close all;

tp    = {'D0' 'D2'  'D5' 'D7' 'D9' 'D12' 'D14'}; %
% Series 1
group1 = {'*A*C*.csv','*A*E*.csv','*B*E*.csv','*B*N*.csv','*B*W*.csv','*F*W*.csv'};
gname1 = {'Control_s1_AC','Control_s1_AE','Control_s1_BE','Control_s1_BN','Control_s1_BW','Control_s1_FW'};

file1 = 'resEnv_sim_comp/K_env/';
%gname = gname1;

%gname = [gname1s, gname1];

% Series 2
%

group2 = {%'*A*C*.csv','*A*N*.csv','*A*S*.csv','*A*W*.csv','*A*E*.csv',...
         %'*B*C*.csv','*B*N*.csv','*B*S*.csv','*B*W*.csv','*B*E*.csv',...
         %'*C*C*.csv','*C*N*.csv','*C*S*.csv','*C*W*.csv','*C*E*.csv',...
         %'*DD*C*.csv','*DD*N*.csv','*DD*S*.csv','*DD*W*.csv','*DD*E*.csv',...
                    '*E*N*.csv','*E*S*.csv','*E*W*.csv'           ,...
                    '*F*N*.csv',           '*F*W*.csv','*F*E*.csv'
         };
%prf_gc2= 'Pac_';
gname2 = {%'Pac_0p5_AC','Pac_0p5_AN','Pac_0p5_AS','Pac_0p5_AW','Pac_0p5_AE',...
          %'Pac_0p05_BC','Pac_0p05_BN','Pac_0p05_BS','Pac_0p05_BW','Pac_0p05_BE',...
          %'Pac_0p005_CC','Pac_0p005_CN','Pac_0p005_CS','Pac_0p005_CW','Pac_0p005_CE',...
          %'Pac_0p0005_DC','Pac_0p0005_DN','Pac_0p0005_DS','Pac_0p0005_DW','Pac_0p0005_DE',...
              'Control_s2_EN','Control_s2_ES','Control_s2_EW'     ,...
              'Control_s2_FN','Control_s2_FW','Control_s2_FE'
          };

%gname2e = append(prf_gc2,gname2);

          
gname = [gname1, gname2]; %, gname1

% File
file2 = 'resEnv_series_2_sim_comp/K_env/';
%


% CA series 1-2
file12s1 = 'resEnv_CA_run14_HYBRID_NO_R_NO_DEATH/'; %'resEnv_CA_s12_run13_HYBRID_NO_R/'; %
%files2   = 'resEnv_CA_run2_TREAT_HYBRID_NO_R/'; 

%group1s = {'*ACs1*.csv','*AEs1*.csv','*BEs1*.csv','*BNs1*.csv','*BWs1*.csv','*FWs1*.csv'};
prf = 'CA_coord_';
%gname1s = {'Control_s1_AC','Control_s1_AE','Control_s1_BE','Control_s1_BN','Control_s1_BW','Control_s1_FW'};
gname1s = append(prf,gname1);
group1s = append('*',gname1,'*.csv');

%file1s = 'resEnv_CA_s1/run4_adhes-5_phenchange0/';

%group2s  = append('*',prf_gc2,gname2,'*.csv');
%group2s = {'*ENs2*.csv','*ESs2*.csv','*EWs2*.csv','*FEs2*.csv','*FNs2*.csv','*FWs2*.csv'};
%gname2s = {'Control_s2_EN','Control_s2_ES','Control_s2_EW','Control_s2_FE','Control_s2_FN','Control_s2_FW'};
gname2s = append(prf,gname2);
group2s = append('*',gname2,'*.csv');
%prf_s2  = 'CA_';
%gname2s = append(prf_s2,prf_gc2,gname2);

%file2s = 'resEnv_CA_s2/K_env/';


gname_ca = [gname1s,gname2s]; %, 

run plotopt.m

figfile = '.'; %'Figures/CA_series_1/run4_adhes-5_phenchange0/';


%%
Kenv = struct;

Kenv = import_Kenv(file1,group1,tp,gname1,Kenv);
disp('---------------')
Kenv = import_Kenv(file2,group2(end-5:end),tp,gname2(end-5:end),Kenv);
disp('---------------')
%Kenv = import_Kenv(file2,group2(1:20),tp,gname2(1:20),Kenv);
disp('---------------')
Kenv = import_Kenv(file12s1,group1s,tp,gname1s,Kenv);
disp('---------------')
Kenv = import_Kenv(file12s1,group2s(end-5:end),tp,gname2s(end-5:end),Kenv);
disp('---------------')
%Kenv = import_Kenv(files2,group2s(1:20),tp,gname2s(1:20),Kenv);
disp('---------------')

%% Plot
c=1;
for i=1:length(group)
    h=figure('Name',gname{i},'Position', [100, 100, 1200, 900],'Visible','Off');
    ah = gobjects(7, 1);
    for j=1:length(tp)
       
        subplot(3,3,c);
        Kenvplot(Kenv.(tp{j}).(gname{i}));
        title(tp{j},'Interpreter','latex');
        c=c+1;
    end
    hh = mtit(' ');
    xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
    set(xlh, 'Visible', 'On');
    ylh=ylabel(hh.ah,'$K(r)-\frac{4}{3}\pi r^3$');
    set(ylh, 'Visible', 'On');
    hL=legend({'Observed','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
    newPosition = [0.4 0.2 0.2 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits);
    c=1;  
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(h,[figfile 'K_' gname{i} '.png'],'-dpng','-r350')
end
%% Envelope from all Paclitaxel 0.5uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p5 = {};
for i=1:length(tp)
   
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.5um
    for j=1:5
        env_0p5{i}{j}=env{i}{j,1};
    end
    subplot(3,3,c);
    envelope(env_0p5{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.5 \mu M$','$SEM_{Pac\ 0.5 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p5.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.05uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p05 = {};
for i=1:length(tp)
    m=1;
    % Paclitaxel 0.05um
    for j=6:10
        env_0p05{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p05{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.05 \mu M$','$SEM_{Pac\ 0.05 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p05.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.005uM  distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p005 = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.005um
    for j=11:15
        env_0p005{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p005{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.005 \mu M$','$SEM_{Pac\ 0.005 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p005.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.0005uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p0005 = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.0005um
    for j=16:20
        env_0p0005{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p0005{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.0005 \mu M$','$SEM_{Pac\ 0.0005 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p0005.png'],'-dpng','-r350')

%% Envelope from all control exp distributions - series 1 and 2
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1100, 900]);
c=1;
env_cont = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=1:12 %21:32
        env_cont{i}{m}=env{i}{j,1};
        m=m+1;
    end
    %length(tp)
    %subplot(3,3,c);
    %envelope(env_cont{i});
    %title(tp{i},'Interpreter','latex');
    %c=c+1;
    
end
%{
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(mm)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'Observed','SEM$_{observed}$','Random','CSR envelope'},'Interpreter','latex','FontSize',24,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.15 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,['K_average_control_all.png'],'-dpng','-r350')
%}

%% Envelope for Pac  exp
env_0p5 = sep_env(env,tp,13,17);
env_0p05 = sep_env(env,tp,18,22);
env_0p005 = sep_env(env,tp,23,27);
env_0p0005 = sep_env(env,tp,28,32);

%% Envelope from control distributions CA series 1 and 2
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1000, 1000]);
c=1;
env_ca = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=13:24 %33:44
        env_ca{i}{m}=env{i}{j,1};
        m=m+1;
    end

    %subplot(3,3,c);
    %envelope(env_cont{i});
    %title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
%%
env_ca_0p5 = sep_env(env,tp,45,49);
env_ca_0p05 = sep_env(env,tp,50,54);
env_ca_0p005 = sep_env(env,tp,55,59);
env_ca_0p0005 = sep_env(env,tp,60,64);


%% Envelope from control distributions series 1 and 2
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1000, 1000]);
c=1;
env_s = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=1:12
        env_s{i}{m}=env{i}{j,1};
        m=m+1;
    end

    %subplot(3,3,c);
    %envelope(env_cont{i});
    %title(tp{i},'Interpreter','latex');
    c=c+1;
    
end

%% Plot K's average for each dataset
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1400, 500]);
h=figure('Name','Average K','Position', [100, 100, 1000,950]);
c=1;

for i=1:length(tp)
    subplot(3,3,c);%,'Position',[0.1+0.99*(i-1)/6 0.2 0.13 0.65]);
    %envelopem({env_ca{i};env_ca_0p0005{i};env_ca_0p005{i};env_ca_0p05{i};env_ca_0p5{i}}); %;env_0p05{i};env_0p005{i};env_0p0005{i}
    %envelopem({env_cont{i};env_0p0005{i};env_0p005{i};env_0p05{i};env_0p5{i}}); 
    envelopem({env_ca{i}; env_cont{i}});
    %envelopem({env_e{end}{1};env_e{end}{2};env_e{end}{3};env_e{end}{4}});
    
    title(tp{i},'Interpreter','latex');
    c=c+1;
    axis tight
end
hh = mtit('','xoff',0.,'yoff',0.09); %'Experiments'
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');%,'FontSize',34
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On'); %,'FontSize',34
%hL=legend({'$Control$','$SEM_{Control}$','$0.0005 \mu M$','$SEM_{0.0005 \mu M}$','$0.005 \mu M$','$SEM_{0.005 \mu M}$','$0.05 \mu M$','$SEM_{0.05 \mu M}$','$0.5 \mu M$','$SEM_{0.5 \mu M}$',...
%    'Random','CSR envelope'},'Interpreter','latex','FontSize',12,...
%        'Location','none','NumColumns',1);
%set(hL,'Position',[0.7 0.15 0.2 0.2]);
hL=legend({'Simulations','SD$_{Sim}$','Experiments','SD$_{Exp}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);     
            %'$SEM_{Control}$',...
            %'$Pac\ 0.5 \mu M$',...
            % '$SEM_{Pac\ 0.5 \mu M}$',...
            %'$Pac\ 0.05 \mu M$',... 
            %'$SEM_{Pac\ 0.05 \mu M}$',...
            %'$Pac\ 0.005 \mu M$',... 
            %'$SEM_{Pac\ 0.005 \mu M}$',...
            %'$Pac\ 0.0005 \mu M$','$SEM_{Pac\ 0.0005 \mu M}$'});
%axis([-inf inf 0 2e9])
%set(hL,'Position',[0.925 0.45 0.04 0.1],'color','none','box','off');
            
newPosition = [0.5 0.2 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,['K_average_all.png'],'-dpng','-r350')
print(h,[file12s1 'K_average_sim_exp_all_cont.png' ],'-dpng','-r350')



%% Test surf
run plotopt.m
h=figure('Name','Average K','Position', [100, 100, 1400, 900]);
%h=figure('Name','Average K','Position', [100, 100, 800,700]);
c=1;
dose={'Control','0.0005 $\mu M$','0.005 $\mu M$','0.05 $\mu M$','0.5 $\mu M$'};

for i=1:length(tp)
    %subplot(2,3,c);
    %envelopem({env_ca{i};env_ca_0p0005{i};env_ca_0p005{i};env_ca_0p05{i};env_ca_0p5{i}}); %;env_0p05{i};env_0p005{i};env_0p0005{i}
    %envelopem({env_cont{i};env_0p0005{i};env_0p005{i};env_0p05{i};env_0p5{i}}); 
    envelope_surf({env_cont{i};env_0p0005{i};env_0p005{i};env_0p05{i};env_0p5{i}},dose);
    %envelopem({env_e{end}{1};env_e{end}{2};env_e{end}{3};env_e{end}{4}});
    
    title(tp{i},'Interpreter','latex');
    c=c+1;
    axis tight
end

newPosition = [0.5 0.2 0.1 0.01];
newUnits = 'normalized';
%set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,['K_average_all.png'],'-dpng','-r350')
%print(h,[files2 'K_average_surf_exp.png' ],'-dpng','-r350')



%% barplots auc control

for i=1:length(tp)

    [au(i,:),su(i,:)]=aucm({env_cont{i},env_ca{i}});

end

cmap = colormap(winter);
%h=figure('Position', [100, 100, 800, 500]);
hold on
hh=bar(au,'EdgeColor','none', 'FaceAlpha',1);
for k = 1:size(hh,2)
    hh(k).FaceColor = cmap(60*k+40,:); %repmat(cmap(k,:),6,1);
end
xdata= get (hh,'XData');
A=xdata{1};
er=errorbar(A-0.15,au(:,1),su(:,1));
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er2=errorbar(A+0.15,au(:,2),su(:,2));
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
hold off

%% barplot p0.5
for i=length(tp):length(tp)

    [au(i,:),su(i,:)]=aucm({env_0p5{i},env_0p05{i},env_0p005{i},env_0p0005{i}});

end
%
cmap = colormap(winter);
%h=figure('Position', [100, 100, 800, 500]);
hold on
hh=bar(au,'EdgeColor','none', 'FaceAlpha',1);
for k = 1:size(hh,2)
    hh(k).FaceColor = cmap(60*k+40,:); %repmat(cmap(k,:),6,1);
end
xdata= get (hh,'XData');
A=xdata{1};
er=errorbar(A-0.15,au(:,1),su(:,1));
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er2=errorbar(A+0.15,au(:,2),su(:,2));
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
hold off



%% Plot treatment experiment and CA together

c=1;
env_e = {};
env_s = {};
for i=1:length(tp)
    m1=1;
    m2=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Experiment
    for j=1:5:20
        c=1;
        for k = j:j+4
        env_e{i}{m1}{c}=env{i}{k,1}; 
        c=c+1;
        end
        m1=m1+1;
    end

    % CA
    for j=21:5:40
        c=1;
        for k = j:j+4
        env_s{i}{m2}{c}=env{i}{k,1}; 
        c=c+1;
        end
        m2=m2+1;
    end
end

%run plotopt.m
set(0,'defaultaxesfontsize',20);
set(0,'defaultaxeslinewidth',1)
set(0,'defaultlinelinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

c=1;
h1=figure('Name','Average K Pac 0.5','Position', [100, 100, 1000, 800]);
for i=1:length(tp)
    subplot(2,3,c);
    envelopem({env_e{i}{1};env_s{i}{1}});
    title(tp{i},'Interpreter','latex');
    c=c+1;
end
hh = mtit('0.5 $\mu m$','yoff',-0.6);
set(hh.th,'edgecolor',.5*[1 1 1]);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On','FontSize',34);
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On','FontSize',34);
hL=legend({'Experiments','SD$_{Exp}$','Simulations','SD$_{Sim}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);            
newPosition = [0.75 0.25 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h1,[file12s 'K_average_sim_vs_exp_Pac0p5.png' ],'-dpng','-r350')

%
c=1;
h2=figure('Name','Average K Pac 0.05','Position', [100, 100, 1000, 800]);
for i=1:length(tp)
    subplot(2,3,c);
    envelopem({env_e{i}{2};env_s{i}{2}});
    title(tp{i},'Interpreter','latex');
     c=c+1;
end
hh = mtit('0.05 $\mu m$','yoff',0.01,'xoff',-0.57);
set(hh.th,'edgecolor',.5*[1 1 1]);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On','FontSize',34);
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On','FontSize',34);
hL=legend({'Experiments','SD$_{Exp}$','Simulations','SD$_{Sim}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);            
newPosition = [0.75 0.25 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h2,[file12s 'K_average_sim_vs_exp_Pac0p05.png' ],'-dpng','-r350')

%

c=1;
h3=figure('Name','Average K Pac 0.005','Position', [100, 100, 1000, 800]);
for i=1:length(tp)
    subplot(2,3,c);
    envelopem({env_e{i}{3};env_s{i}{3}});
    title(tp{i},'Interpreter','latex');
     c=c+1;
end
hh = mtit('0.005 $\mu m$','yoff',0.02,'xoff',-0.58);
set(hh.th,'edgecolor',.5*[1 1 1]);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On','FontSize',34);
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On','FontSize',34);
hL=legend({'Experiments','SD$_{Exp}$','Simulations','SD$_{Sim}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);            
newPosition = [0.75 0.25 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h3,[file12s 'K_average_sim_vs_exp_Pac0p005.png' ],'-dpng','-r350')

%
c=1;
h4=figure('Name','Average K Pac 0.0005','Position', [100, 100, 1000, 800]);
for i=1:length(tp)
    subplot(2,3,c);
    envelopem({env_e{i}{4};env_s{i}{4}});
    title(tp{i},'Interpreter','latex');
     c=c+1;
end
hh = mtit('0.0005 $\mu m$','yoff',0.02,'xoff',-0.57);
set(hh.th,'edgecolor',.5*[1 1 1]);
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On','FontSize',34);
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On','FontSize',34);
hL=legend({'Experiments','SD$_{Exp}$','Simulations','SD$_{Sim}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);            
newPosition = [0.75 0.25 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h4,'Units','Inches');
pos = get(h4,'Position');
set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h4,[file12s 'K_average_sim_vs_exp_Pac0p0005.png' ],'-dpng','-r350')




%% Calculate differences of K between two subsequent days
%{
for i=1:length(tp)-1
    j=i+1;
    for k=1:length(group)
   
        dKenv.([tp{i} tp{j}]).(gname{k})=[Kenv.(tp{j}).(gname{k}).r, Kenv.(tp{j}).(gname{k}).obs - Kenv.(tp{i}).(gname{k}).obs];
        
    end
end

%% Plot the envelope of the differences
h=figure('Name','Average dK','Position', [100, 100, 1200, 900]);
for i=1:length(tp)-1
   
    j=i+1;
    denv{i}=struct2cell(dKenv.([tp{i} tp{j}]));
    subplot(3,3,i);
    diffenvelope(denv{i});
    title([tp{j} '-' tp{i}],'Interpreter','latex');
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle \Delta K(r) \rangle $');
set(ylh, 'Visible', 'On');
hL=legend({'Observed','Observed SEM'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%        'Location','southeastoutside','NumColumns',1);
newPosition = [0.8 0.2 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,['Figures/diffK_average.png'],'-dpng','-r0')
%}
