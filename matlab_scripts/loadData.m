clear data
clear
%change the directory here to load other samples, like if you want to load
%slide1, make relPath = 'slide1';
relPath = pwd; 
files = [dir('../Dan_data/vesicle_combined/x*.xml')];
for i = 1:length(files)
   data(i).name = files(i).name(1:end-4); 
   [data(i).mz, data(i).intens,~] = parseSolarixXML(fullfile(files(i).folder, files(i).name));
end
%%
%and make this sample specific
save '../Dan_data/vesicle_combined_aligned.mat' data
%save '../data/Mito_all_files.mat' files
%%
mzVals = [153 242 345 442 504 603 747 834 906 1061 1544];
bins = [.0004 .0007 .001 .002 .004 .006 .007 .008 .01 .02 .03]/2;
mzLimits = [150, 1500];
data = alignMasses({'../Dan_data/vesicle_combined.mat'},mzVals,bins,mzLimits);