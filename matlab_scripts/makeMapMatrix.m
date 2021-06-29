%create a map matrix to check the compound coverage of each cell with
%singals detected

for i=1:map_dim(1)
    cell_id = cat(2,cell_id,match_map{i,7});
end
cell_id = unique(sort(cell_id,'ascend'));

map_matrix = zeros(length(category),length(cell_id));

for i=1:map_dim(1)
    cell_list = match_map{i,7};
    [r1,c1] = find(strcmp(category,match_map(i,1)));
    for j=1:length(cell_list)
        [r,c] = find(cell_id == cell_list(1,j));
        
        map_matrix(c1,c) = 1;
    end
end

figure();
imagesc(map_matrix)
%set(gca, 'XTick', 1:length(cell_id)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(category)); % center y-axis ticks on bins
%set(gca, 'XTickLabel', cell_id); % set x-axis labels
set(gca,'YTickLabel', category); % set y-axis labels
%colormap('jet');