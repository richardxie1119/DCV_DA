function [sum_int,mz] = sum_spectra(data)

dim = length(data);

sum_int=zeros(1,dim);
for i=1:dim
    
    sum_int = sum_int + data(i).intens;
    mz = data(i).mz;
end


end