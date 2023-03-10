function plot_centr(group1,gname1,coord1,tp,name)

for i=1:length(group1)
    
%    f=figure('Name',gname1{i},'Position', [100, 100, 1400, 1200]);
    NameArray = {'MarkerSize'};
    ValueArray = {3};
    for j=1:length(tp)
        A=coord1.(gname1{i}).(tp{j});
        f=figure('Name',gname1{i},'Position', [100, 100, 800, 600]);
%        h(j)=subplot(3,3,j);
%        hax=gca;
%        hax.FontSize=16;
%        cmap = jet(length(A(:,3)));
         AxesH = axes('Parent', f, ...
             'NextPlot', 'add');  % Equivalent to: "hold on"
         xlabel('x $(\mu m)$','FontSize',32);
         ylabel('y $(\mu m)$','FontSize',32);
         zlabel('z $(\mu m)$','FontSize',32);
         title(tp{j},'Interpreter','latex','FontSize',36);
         axis([0 2500 0 2500 0 917])

%        for k=1:length(A(:,3))
%            plot3(A(k,1),A(k,2),A(k,3),'o','color',cmap(k,:),'markerfacecolor',cmap(k,:))
%        end
        scatter3(A(:,1),A(:,2),A(:,3),25,A(:,3),'filled');
        view(50,35)
        colormap(jet);
%        A(j)=copyobj(allchild(get(A(j),'CurrentAxes')),h(j));
%        view(20,20)
        %view(2)
%        set(A(j),NameArray,ValueArray)
      
        %hold off
        clear A
        set(f,'Units','Inches');
        pos = get(f,'Position');
        set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        disp('Saving...')
        print(f,[name '_' gname1{i} '_tp_' num2str(tp{j}) '.png'],'-r350','-dpng')
    end
    %subplot(3,3,8);
    %hcx=gca;
    %hcx.FontSize=16;
    %hold on
    %plot(time,struct2array(count.(gname{i})),'.-','LineWidth',4,'MarkerSize',30);
    %errorbar(time,gcm,gcstd,'.-','LineWidth',4,'MarkerSize',30)
    %xlabel('time $(day)$','Interpreter','latex')
    %ylabel('Nuclei count','Interpreter','latex')
    %hold off
    
%    set(f,'Units','Inches');
%    pos = get(f,'Position');
%    set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%    disp('Saving...')
%    print(f,['centroids_' gname1{i} '.png'],'-r350','-dpng')
end


end