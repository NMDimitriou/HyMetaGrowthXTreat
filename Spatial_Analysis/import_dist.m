function [IND,KND] = import_dist(path,name1,name2,group1,gname1,tp,IND,KND)


for i=1:length(group1)
    
    disp(['Opening ' group1{i}])
    sampIND{i}=dir([path name1 group1{i}]);
    sampKND{i}=dir([path name2 group1{i}]);

    namesIND     = {sampIND{i}.name};
    namesIND     = natsort(namesIND);
    namesKND     = {sampKND{i}.name};
    namesKND     = natsort(namesKND);

    for j=1:length(tp)

        k=find(contains(namesIND,tp{j}));
        disp(['   ->' namesIND{k}])
        IND.(gname1{i}).(tp{j})=readmatrix([path namesIND{k}]);
        KND.(gname1{i}).(tp{j})=readmatrix([path namesKND{k}]);   
    end
    
end


end