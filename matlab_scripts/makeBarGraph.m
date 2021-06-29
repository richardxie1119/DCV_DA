%plot the bar graph
map_dim = size(match_map);
category = {};
int_vals = [];
err_vals = [];
cell_num = [];
cell_id = [];
index=1;

for i=1:map_dim(1)
    temp = char(match_map{i,1});
    if ~any(strcmp(category,temp))
        category{1,index} = char(match_map{i,1});
        int_vals(1,index) = match_map{i,2};
        err_vals(1,index) = match_map{i,5};
        cell_num(1,index) = match_map{i,6};
        index = index+1;
    
    else
        [r,c] = find(strcmp(category, temp));
        %category{1,i} = strcat(char(match_map{i,1}),'_replicate');
        int_vals(1,c) = int_vals(1,c)+match_map{i,2};
        err_vals(1,c) = err_vals(1,c)+match_map{i,5};
        cell_num(1,c) = cell_num(1,c)+match_map{i,6};
    end


end

species = categorical(category);
ppm_vals = err_vals.*10^6;

figure();
subplot(3,1,1);bar(species,int_vals)
title('sum intensity values');
subplot(3,1,2);bar(species,ppm_vals)
title('ppm of summed signals')
subplot(3,1,3);bar(species,cell_num)
title('total cell number')


