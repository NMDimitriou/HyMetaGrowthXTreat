function varargout = envelope(A)
%% Plots the envelope of the summary function for a three-dimensional point pattern
% from all the samples.
% input: the matrix obtained from envelope.pp3 function of spatstat package
% for all samples
% of R
% output: the plot
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

x=A{1}.r';
ytheo=A{1}.theo'-4*pi*x.^3/3;

for i=1:length(A)
    
    yobs_all(:,i) = A{i}.obs;
    ylo_all(:,i)  = A{i}.lo;
    yhi_all(:,i)  = A{i}.hi;
    
end

%for i=1:length(yobs_all(:,1))
   
    %yobs_min(i)=min(yobs_all(i,:))'-4*pi*x(i).^3/3;
    %yobs_max(i)=max(yobs_all(i,:))'-4*pi*x(i).^3/3;
    
    %ylo_min(i)=min(ylo_all(i,:))'-4*pi*x(i).^3/3;
    %yhi_max(i)=max(yhi_all(i,:))'-4*pi*x(i).^3/3;
    
%end

yobs_mean = mean(yobs_all,2)' -4*pi*x.^3/3;
yobs_std  = std(yobs_all,0,2)'./sqrt(length(yobs_all(1,:)));
ylo_min  = min(ylo_all,[],2)'  -4*pi*x.^3/3;
yhi_max  = max(yhi_all,[],2)'  -4*pi*x.^3/3;

ycsr_env  = [ylo_min;   yhi_max];
%yobs_env  = [yobs_min; yobs_max];
yobs_ms   = [yobs_mean-yobs_std; yobs_mean+yobs_std];


%figure('Name',name); 
hold on
hax=gca;
plot      (x./1000,yobs_mean,'--','LineWidth',4)
plotshaded(x./1000,yobs_ms ,'b')
plot      (x./1000,ytheo,'--','LineWidth',3)
plotshaded(x./1000,ycsr_env ,[17 17 17]/255)
%plotshaded(x,yobs_env ,'b')
%xlabel('Neighborhood radius $r$ $(\mu m)$','Interpreter','latex')
%ylabel('$K(r)-\frac{4}{3}\pi r^3$','Interpreter','latex')
%legend({'Observed','Random','CSR envelope'},'Interpreter','latex','FontSize',14,...
%         'Location','northwest','NumColumns',1)
hax.FontSize=24;
axis([0 2 -1.5e+09 1.5e+09]) %max(x./1000)
box off