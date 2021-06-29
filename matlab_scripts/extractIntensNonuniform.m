function [intens] = extractIntensNonuniform(mz, data, bins, mzRange)

    function [ denom ] = extractIntensHelp( mzList, scan, bin )
        idx = rangesearch(scan.mz', mzList, bin/2);
        denom = cellfun(@(a) sum(scan.intens(a)), idx);
    end
    mz = mz';
    intens = [];
    for i = 1:length(bins)
        t = mz > mzRange(i,1) & mz < mzRange(i,2);
        t = arrayfun(@(a) extractIntensHelp(mz(t), a, bins(i)), data(:),...
            'uniformoutput', false);
        t = [t{:}];
        intens = [intens; t];
    end
end