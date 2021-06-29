function [numerator, denomenator,spec] = centOfMassNonuniform(mz, data, bins, mzRange)

    numerator = [];
    denomenator = [];
    
    %concatenate all spectra together
    spec = arrayfun(@(a) a.mz, data(:), 'uniformoutput', false);
    intens = arrayfun(@(a) a.intens, data(:), 'uniformoutput', false);
    spec = [horzcat(spec{:}); horzcat(intens{:})]';
    mz = mz';
    
    for i = 1:length(bins)
        t = mz > mzRange(i,1) & mz < mzRange(i,2);
        %find idx
        idx = rangesearch(spec(:,1), mz(t), bins(i)/2);
        %calculate sum and sumsquare
        n = cellfun(@(a) sum(prod(spec(a,:),2)), idx);
        d = cellfun(@(a) sum(spec(a,2)), idx);
        numerator = [numerator; n];
        denomenator = [denomenator; d];
    end
end

