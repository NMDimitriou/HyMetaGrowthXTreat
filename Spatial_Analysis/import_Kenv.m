function Kenv = import_Kenv(file2,group2,tp,gname2,Kenv)

for j=1:length(tp)
    %samp2=dir([file2 tp{j}]); % {j}
    %names   = {samp2{i}.name};
    for i=1:length(group2)
        samp2{i,j}=dir([file2 tp{j} '/' group2{i}]); % {j}
        disp(['Opening '  samp2{i,j}.name ' -> ' gname2{i} ', ' tp{j}])
        Kenv.(tp{j}).(gname2{i})=readtable([samp2{i,j}.folder '/' samp2{i,j}.name]);
    end
end

end