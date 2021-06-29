function [match_map] = matchMassList(data,mass_list_dir,mode,precision)

[num,txt,raw] = xlsread(mass_list_dir,mode);
num(:,1) = [];
raw(:,2) = [];
match_map={};
index = 1;
dim = size(num);

for i=1:length(data.mzs)
    int_sum = sum(data.intens(i,:));
    total_cell = nnz(data.intens(i,:));
    cell_id = find(data.intens(i,:) ~=0);
    if int_sum ~= 0
        mz = data.mzs(i);
        [r,c] = find(num-rem(num,precision) == mz-rem(mz,precision));
        match_temp = [r,c];
        match_dim = size(match_temp);
        %disp(length(r))
        %disp(match_temp)
        if isempty(match_temp)
            %
        else
            for j=1:match_dim(1)
                match_map{index,1} = char(raw(match_temp(j,1)+1,1));
                match_map{index,2} = int_sum;
                match_map{index,3} = mz;
                match_map{index,4} = num(match_temp(j,1),match_temp(j,2));
                match_map{index,5} = abs(mz-num(match_temp(j,1),match_temp(j,2)))/num(match_temp(j,1),match_temp(j,2));
                match_map{index,6} = total_cell;
                match_map{index,7} = cell_id;
                match_map{index,8} = char(raw(1,c+1));
                index=index+1;
            end
           
        end
    end
end

%plot the bar graph
% map_dim = size(match_map);
% category = {};
% int_vals = [];
% for j=1:map_dim(1)
%     category{1,j} = match_map{j,1};
%     int_vals(1,j) = match_map{j,2};
% end
% species = categorical(category);
% 
% figure();
% bar(species,int_vals)
