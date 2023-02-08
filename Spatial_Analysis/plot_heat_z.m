function plot_heat_z(coord,ca_coord1,tp,zmax,n,name,tit)

gc=struct2array(coord);
gc=struct2table(gc);
gc=table2array(gc);
gc_ca=struct2array(ca_coord1);
gc_ca=struct2table(gc_ca);
gc_ca=table2array(gc_ca);
[~,cols]=size(gc);
%[l2,d2]=size(gc_ca);
%h=15;


z = linspace(0,zmax,n);

heat_sim = zeros(cols,n-1);
heat_exp = zeros(cols,n-1);

for i=1:cols
    
    tmp=vertcat(gc{:,i});
    tmp=tmp(:,3);
    tmp_ca=vertcat(gc_ca{:,i});
    tmp_ca=tmp_ca(:,3);
    
    for j=1:n
        
        
        [heat_exp(i,:),~] = histcounts(tmp,z);
        [heat_sim(i,:),~] = histcounts(tmp_ca,z);
        
        heat_exp(i,:) = heat_exp(i,:)./sum(heat_exp(i,:));
        heat_sim(i,:) = heat_sim(i,:)./sum(heat_sim(i,:));
    end
end

[T,Z]=meshgrid(tp,z(1:end-1));



f1=figure('Position', [100, 100, 500, 500],'Visible','on');
%ax(1)=subplot(1,2,1);
surf(T,Z,heat_exp','FaceAlpha',0.7,'EdgeColor','none') %,'EdgeColor','none'
colormap(redblue)
view(2)
title(['Experiments']) % - ' tit
axis tight
lims = clim;
h=colorbar ;
xlabel('time (days)');
ylabel('$z\ (\mu m)$');
annotation('textbox',...
     h.Position,...
    'FitBoxToText','off',...
    'FaceAlpha',0.3,...
    'BackgroundColor',[1 1 1]);

grid off


f2=figure('Position', [100, 100, 500, 500],'Visible','on');
%ax(2)=subplot(1,2,2);
surf(T,Z,heat_sim','FaceAlpha',0.7,'EdgeColor','none') %,'EdgeColor','none'
colormap(redblue)
view(2)
title(['Simulations']) % - ' tit
axis tight
caxis(lims);
h=colorbar ;
xlabel('time (days)');
ylabel('$z\ (\mu m)$');
annotation('textbox',...
     h.Position,...
    'FitBoxToText','off',...
    'FaceAlpha',0.3,...
    'BackgroundColor',[1 1 1]);

grid off

%{
h=colorbar;
set(h, 'Position', [.8714 .181 .0281 .6150])
for i=1:2
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [pos(1) pos(2) 0.8*pos(3) pos(4)]);
end
annotation('textbox',...
     h.Position,...
    'FitBoxToText','off',...
    'FaceAlpha',0.3,...
    'BackgroundColor',[1 1 1]);

hh = mtit(' ');
xlh=xlabel(hh.ah,'time (days)');
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$z\ (\mu m)$');
set(ylh, 'Visible', 'On')
colormap(redblue)
%}
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,[name  '_1.png'],'-r350','-dpng')

set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f2,[name  '_2.png'],'-r350','-dpng')


end