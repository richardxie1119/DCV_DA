[num,txt,raw] = xlsread('../mass_calculation','Negative');
num(:,1) = [];
raw(:,2) = [];
compounds_list = raw(2:end,1);
type = raw(1,3:end);
match_map={};
index = 1;
dim = size(num);

data_size = size(data);

for i=1:data_size(2)
    mzs = data(i).mz;
    intens = data(i).intens;
    match_compounds = {};
    match_ppm = [];
    match_type = {};
    for j=1:length(type)
        mass_list = num(:,j);
        [a,b]=ismembertol(mzs,mass_list,0.005,'DataScale',1);
        if sum(a)== 0
            %
        else 
            mzs_measured = mzs(a);
            mass_list_index = b(a);
            for k=1:length(mzs_measured)
                match_compounds{1,k} = compounds_list(mass_list_index(k));
                match_ppm(k) = (mzs_measured(k)-mass_list(mass_list_index(k)))/mass_list(mass_list_index(k));
                match_type{1,k} = type(j);
            end
        end
    end
    
    match_map{i,1} = match_compounds;
    match_map{i,2} = match_ppm;
    match_map{i,3} = match_type;
    
end

                