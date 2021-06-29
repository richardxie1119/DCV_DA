mzVals = [150,200,250,300,350,400,450,500, 550, 600, 650, 700, 750, 800,850,900,950,1100];
bins = [0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035, 0.004, 0.0045, 0.005, 0.0055, 0.006, 0.0065,0.007,0.0075,0.008,0.0085,0.009];
mzLimits = [150, 1100];

%get initial mz list
binDiv = 18;
fineBins = mzLimits(1):bins(1)/binDiv:(mzVals(1) + mzVals(2))/2;
mzRange = [mzLimits(1) (mzVals(1) + mzVals(2))/2];
for i = 2:length(mzVals)-1
    fineBins = [fineBins (mzVals(i-1) + mzVals(i))/2:bins(i+1)/binDiv:(mzVals(i) + mzVals(i+1))/2];
    mzRange = [mzRange ; (mzVals(i-1) + mzVals(i))/2 (mzVals(i) + mzVals(i+1))/2];
end
fineBins = [fineBins (mzVals(end-1) + mzVals(end))/2:bins(i+1)/binDiv:mzLimits(2)];
fineBins = unique(fineBins);
mzRange = [mzRange ; (mzVals(end-1) + mzVals(end))/2 mzLimits(2)];

matFiles_names = dir('*.mat');
matFiles = dir('*.mat');
%remove mat files generated from previous runs, may be overwritten!
matFiles = matFiles(cellfun(@isempty, strfind({matFiles(:).name}, '_nonU')));

%need to find counts of all images separately
counts1 = zeros(1,length(fineBins)-1);
cells = 0;
disp('Populating m/z histogram')
for i = 1:length(matFiles)
    disp(matFiles_names(i).name);
    file_name = matFiles_names(i).name;
    load(file_name)
    %propagate master mass list
    mzs = arrayfun(@(a) transpose(a.mz), data, 'uniformoutput', false);
    mzs = vertcat(mzs{:});
    counts = counts1 + histcounts(mzs, fineBins);
    cells = cells + length(data);
end

fineBins = (fineBins(1:end-1) + fineBins(2:end))/2; %move bin to center of range

%find peaks in "histogram", exclude those in less than 0.1% of cells
figure; plot(fineBins, counts);
hold on; plot([fineBins(1) fineBins(end)], [cells*0.001 cells*0.001]);
mzList = [];
%this will not handle peaks at the boundaries of the ranges
for i = 1:length(bins)
    t = fineBins > mzRange(i,1) & fineBins < mzRange(i,2);
    [~, t] = findpeaks(counts(t), fineBins(t),...
        'MinPeakHeight', cells*0.001,...
        'MinPeakDistance', bins(i));
    mzList = [mzList t];
end


%reload each image to get their contribution to the mass list
numerator = zeros(size(mzList));
denomenator = zeros(size(mzList));
disp('Calculating center of mass')
for i = 1:length(matFiles)
    disp(matFiles_names(i).name);
    file_name = matFiles_names(i).name;
    load(file_name)
    %propagate master mass list
    [n,d,spec] = centOfMassNonuniform(mzList, data, bins, mzRange);
    numerator = numerator + n';
    denomenator = denomenator + d';
end

%calculate center of mass
mzs = numerator ./ denomenator;
disp(max(abs(mzs - mzList)))
% mzList = mzs;

%record intensity matrix, mzList, save file
disp('Recalculating intensity matrix')
clear d
for j = 1:length(matFiles)
    disp(matFiles_names(j).name);
    file_name = matFiles_names(j).name;
    load(file_name)
    %record new peak picked intensities
    d.names = {data(:).name};
    d.mzs = mzList;
    d.intens = extractIntensNonuniform(mzList, data, bins, mzRange);
    data = d;
    %save
    %save([ matFiles(j).name(1:end-4) '.mat'], 'data')
end
save 'out/pk_data_aligned_800_05.mat' data
