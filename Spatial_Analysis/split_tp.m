function ca_coord1_tsp = split_tp(groupc1, gname1, ca_coord1, ts, iter)

for i=1:length(groupc1)
    %for k=1:iter
        B = [];
        c=1;
        len=length(ca_coord1.(gname1{i})(:,1));%{k}
        for j=2:len

            if(ca_coord1.(gname1{i})(j,1)==ca_coord1.(gname1{i})(j-1,1))%{k}
                B=[B;ca_coord1.(gname1{i})(j-1,:)];%{k}
                ca_coord1_tsp.(gname1{i}).(ts{c})=B; %.(['iter_' num2str(k)])
            else
                c=c+1;
                disp(['c = ' num2str(c), ', len = ' num2str(len) ', j = ' num2str(j)])
                disp(['Changing to ts = ' ts{c}]);
                B=[];
                ca_coord1_tsp.(gname1{i}).(ts{c})=B; %.(['iter_' num2str(k)])
            end
        end
    %end
end

end