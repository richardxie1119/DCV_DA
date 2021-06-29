function [P PFWHH PEXT] = mspeaks_(X,Y,varargin)
%MSPEAKS robust peak detection in noisy signals
%
%  P = MSPEAKS(X,Y) finds relevant peaks in a noisy signal by following
%  three steps: 1) the signal is smoothed using the undecimated wavelet
%  transform with the Daubechies filter banks, 2) putative peak locations
%  and their extents are found, both are resolved beyond the resolution of
%  the scale in X, and 3) apply multiple post-filtering criteria that
%  reduces oversegmented and noisy peaks.
%
%  X and Y are column or row vectors where paired values represent points
%  in the signal. Y can be a matrix with several signals, either
%  column-wise or row-wise, all sharing the same X scale. Units in the X
%  scale (separation units or s.u.) may quantify wavelength, frequency,
%  distance, time or m/z depending on the type of instrument that generates
%  the signal.
%
%  The output peak list P is a matrix with two columns with peak
%  information, the first column has the location of the peaks in the
%  separation axis and the second column has the intensity of the peaks.
%  When more than one signal is given in Y, P is a cell array containing
%  multiple peak lists.
%
%  [P,PFWHH] = MSPEAKS(X,Y) returns a matrix with two columns indicating
%  the left and right location of the full width at half height markers for
%  every peak. For any peak not resolved at FWHH, the peak shape extents
%  are returned instead. PFWHH is a cell array of matrices when several
%  signals are given in Y.  
%
%  [P,PFWHH,PEXT] = MSPEAKS(X,Y) returns a matrix with two columns
%  indicating the left and right location of the peak shape extents
%  determined after wavelet denoising. PEXT is a cell array of matrices
%  when several signals are given in Y. The peak extents are either the
%  closest valley minima or the zero crossing to each side of the peak.
%
%  MSPEAKS(...,'BASE',WB) selects the wavelet base. WB is an integer
%  between 2 and 20. WB defaults to 4.
%
%  MSPEAKS(...,'LEVELS',WL) selects the number of levels for the wavelet
%  decomposition. WL is an integer between 1 and 12. WL defaults to 10.
%
%  MSPEAKS(...,'NOISEESTIMATOR',NE) method to estimate the threshold (T)
%  to filter out noisy components in the first high band decomposition.
%  Options are: 
%
%    'mad'       - Median absolute deviation: 
%    (default)     T = sqrt(2*log(n))*mad(y_h)/.6745
%
%    'std'       - Standard deviation: T = std(y_h)
%
%    T           - User specified threshold, T is a positive real scalar.
%
%  MSPEAKS(...,'MULTIPLIER',M) selects a threshold multiplier constant. M
%  is a positive real scalar. By default M = 1.
%
%  MSPEAKS(...,'DENOISING',FALSE) disables wavelet denoising. Use this
%  option to find peaks from signals that are already smoothed, e.g. with
%  MSLOWESS or MSSGOLAY. Default is true.
%
%  MSPEAKS(...,'PEAKLOCATION',R) sets the ratio of the peak height that
%  selects the points used to compute the centroid mass of the respective
%  peak. R is a scalar value between 0 and 1. Observe that when R=1 the
%  peak location occurs exactly at the maximum of the peak, while when R=0
%  the peak location is computed with all the points from the closest
%  minimum to the left of the peak to the closest minimum to the right of
%  the peak. R defaults to 1.
%
%  MSPEAKS(...,'FWHHFILTER',FF) sets the minimum permissible full width at
%  half height for reported peaks, all others are removed from the output
%  peak list. FF is a positive real scalar given in s.u. and defaults to 0.
%  
%  MSPEAKS(...,'OVERSEGMENTATIONFILTER',OF) sets the minimum permissible
%  distance between neighbor peaks. When the signal is not appropriately
%  smoothed, multiple maxima can appear representing the same peak. By
%  setting OF, oversegmented peaks are joined into a single one. OF is a
%  positive real scalar given in s.u. and defaults to 0.
%
%  MSPEAKS(...,'HEIGHTFILTER',HF) sets the minimum permissible
%  height for reported peaks, all others are removed from the output
%  peak list. FF is a positive real scalar and defaults to 0.
%
%  MSPEAKS(...,'SHOWPLOT',SP) plots the original and the smoothed signal
%  with a mark at the output peak locations. When SP is TRUE, the first
%  signal in Y is used. If MSPEAKS is called without output arguments, a
%  plot will be shown unless SP is FALSE. SP can also contain an index to 
%  one of the signals in Y.
%
%  MSPEAKS(...,'STYLE',S) selects the style for marking the peaks on the
%  plot. Options are:
%     'peak'(default) - places a marker at the crest
%     'exttriangle'   - draws a triangle using the crest and the extents
%     'fwhhtriangle'  - draws a triangle using the crest and the FWHH points
%     'extline'       - places a marker at the crest and vertical lines at
%                       the extents 
%     'fwhhline'      - places a marker at the crest and a horizontal line
%                       at FWHH 
%
%  Examples:
% 
%      load sample_lo_res
%
%      % Adjust the baseline of the eight spectra stored in Y_lo_res:
%      YB = msbackadj(MZ_lo_res,Y_lo_res);
% 
%      % Convert the raw mass spectrometry data to a peak list by finding
%      % the relevant peaks in each spectrum:
%      P = mspeaks(MZ_lo_res,YB);
%
%      % Plot the third spectrogram in YB, the baseline-corrected
%      % intensity values, with the detected peaks marked. 
%      P = mspeaks(MZ_lo_res,YB,'SHOWPLOT',3);
%
%      % Smooth the signal using MSLOWESS (or MSSGOLAY) and then find the
%      % peaks with MSPEAKS, now Wavelet denoising is not necessary:
%      YS = mslowess(MZ_lo_res,YB,'SHOWPLOT',3);
%      P = mspeaks(MZ_lo_res,YS,'DENOISING',false,'SHOWPLOT',3);
%
%      % Use the CELLFUN function to remove all peaks with m/z values less
%      % than 2000 from the eight peaks lists in output P. Then plot the
%      % peaks of the third spectrum (in red) over its smoothed signal (in
%      % blue):   
%      Q = cellfun(@(p) p(p(:,1)>2000,:),P,'UniformOutput',false);
%      figure
%      plot(MZ_lo_res,YS(:,3),'b',Q{3}(:,1),Q{3}(:,2),'rx')
%      xlabel('Mass/Charge (M/Z)')
%      ylabel('Relative Intensity')
%      axis([0 20000 -5 95])
% 
%  See also LCMSDEMO, MSBACKADJ, MSDOTPLOT, MSLOWESS, MSPALIGN,
%  MSPREPRODEMO, MSPPRESAMPLE, MSSGOLAY. 
 
%   Copyright 2006-2012 The MathWorks, Inc.

 
% References:
% [1] J.S. Morris, et.al. "Feature extraction and quantification for mass
%     spectrometry in biomedical applications using the mean spectrum",
%     Bioinformatics 21(9):1764, 2005 
%
% [2] Y. Yasui, et.al. "A data-analytic strategy for protein biomarker
%     discovery: profiling of high-dimensional proteomic data for cancer
%     detection", Biostatistics 4:449-463, 2003
%
% [3] D. Donoho and I. Johnstone, "Adapting to unknown smoothness via
%     wavelet shrinkage," J. Am. Statist. Asso. 90:1200-1224, 1995.
%
% [4] Gilbert Strang and Truong Nguyen, "Wavelets and Filter Banks"
%     Wellesley Cambridge Press, 1996. 
%
% [5] K.R. Coombes, et.al. "Improved peak detection and quantification of
%     mass spectrometry data acquired from surface-enhanced laser
%     desorption and ionization by denoising spectra with the undecimated
%     discrete wavelet transform", Proteomics 5(16):4107-4117, 2005.


% check inputs
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

bioinfochecknargin(nargin,2,mfilename);
% set defaults
thresholdMultiplier = 3;
noiseEstimator = 'mad';
noiseEstimatorConstant = 1;
waveletLevels = 10;              
waveletBasis = 4;
waveletEnable = true;
peakLocationRatio = 1;
fwhhFilter = 0;
oversegFilter = 0;
heightFilter = 0;
style = 'peak';
if nargout == 0
    plotId = 1; 
else
    plotId = 0;
end

if  nargin > 2
    if rem(nargin,2) == 1
        error(message('bioinfo:mspeaks:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'base','levels','noiseestimator','multiplier','denoising',...
              'peaklocation','fwhhfilter','oversegmentationfilter',...
              'heightfilter','showplot','style'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mspeaks:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mspeaks:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  %'base'
                    if ~isnumeric(pval) || ~isscalar(pval) || rem(pval,1) || pval<2 || pval>20
                        error(message('bioinfo:mspeaks:InvalidBase'));
                    else
                        waveletBasis = pval;
                    end
                case 2  %'levels'
                    if ~isnumeric(pval) || ~isscalar(pval) || rem(pval,1) || pval<1 || pval>12
                        error(message('bioinfo:mspeaks:InvalidLevels'));
                    else
                        waveletLevels = pval;
                    end
                case 3  %'noiseestimator'
                    if isnumeric(pval) && isscalar(pval) && pval>0
                        noiseEstimatorConstant = pval;
                        noiseEstimator = 'none';
                    else
                        noiseEstimators = {'mad','std','none'};
                        noiseEstimator = strmatch(lower(pval),noiseEstimators); 
                        if isempty(noiseEstimator) 
                            error(message('bioinfo:mspeaks:NotValidNoiseEstimator'))
                        end
                        noiseEstimator = noiseEstimators{noiseEstimator};
                    end
                case 4  %'multiplier'
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<=0
                        error(message('bioinfo:mspeaks:InvalidMultiplier'));
                    else
                        thresholdMultiplier = pval;
                    end
                case 5  %'denoising'
                    waveletEnable = bioinfoprivate.opttf(pval,okargs(k),mfilename);
                case 6  %'peaklocation'
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0 || pval>1
                        error(message('bioinfo:mspeaks:InvalidPeakLocation'));
                    else
                        peakLocationRatio = pval;
                    end
                case 7  %'fwhhfilter'
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0
                        error(message('bioinfo:mspeaks:InvalidFwhhFilter'));
                    else
                        fwhhFilter = pval;
                    end
                case 8  %'oversegmentationfilter'
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0
                        error(message('bioinfo:mspeaks:InvalidOversegFilter'));
                    else
                        oversegFilter = pval;
                    end
                case 9  %'heightfilter'
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0
                        error(message('bioinfo:mspeaks:InvalidHeightFilter'));
                    else
                        heightFilter = pval;
                    end
                case 10 %'showplot'
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:mspeaks:SPNoScalar'))
                            end 
                        else
                            plotId = 1;
                        end
                    else
                        plotId = 0;
                    end
                case 11 %'style'
                    styles = {'peak','exttriangle','fwhhtriangle','extline','fwhhline'};
                    style = strmatch(lower(pval),styles); 
                    if isempty(style) 
                        error(message('bioinfo:mspeaks:NotValidStyle'))
                    elseif length(style)>1
                        error(message('bioinfo:mspeaks:AmbiguousStyle'))
                    end
                    style = styles{style};
            end
        end
    end
end

% validate X and Y

if ~isnumeric(Y) || ~isreal(Y) || ndims(Y)>2
   error(message('bioinfo:mspeaks:IntensityNotNumericAndReal')) 
end

n = numel(X);

if ~isnumeric(X) || ~isreal(X) || ~isvector(X) || any(diff(X)<=0) || n<2
   error(message('bioinfo:mspeaks:XNotNumericAndReal')) 
end

if diff(size(X))>0
    X = X(:);
    switch n
        case size(Y,2)
            Y = Y';
        case size(Y,1)
            warning(message('bioinfo:mspeaks:XYrowLengthIncompatible'))
    end
else
    switch n
        case size(Y,1)
        case size(Y,2)
            Y = Y';
            warning(message('bioinfo:mspeaks:XYcolumnLengthIncompatible'))
    end
end

[nY,numSignals] = size(Y);

if nY ~= n
    error(message('bioinfo:mspeaks:NotEqualNumberOfSamples'))    
end

P = cell(numSignals,1);
if nargout>1
    PFWHH = cell(numSignals,1);
end
if nargout>2
    PEXT = cell(numSignals,1);
end

    

if plotId>0 %save y for later plotting when necessary
    y_n = Y(:,plotId);
end

% perform wavelet decomposition
if waveletEnable
    if nargout==0 && plotId>0
        Y(:,plotId) = waveletdenoise(Y(:,plotId),waveletBasis,waveletLevels,noiseEstimator,thresholdMultiplier.*noiseEstimatorConstant);
    else
        Y = waveletdenoise(Y,waveletBasis,waveletLevels,noiseEstimator,thresholdMultiplier.*noiseEstimatorConstant);
    end
end

for i = 1:numSignals
if nargout>0 || (i == plotId) 

    y = Y(:,i);
    
    % robust valley finding
    h = [0;find(diff(y));n];
    g = diff(y(h([2 2:end])))<=0 & diff(y(h([2:end end])))>=0;
    leftMin = h([g;false])+1;
    rightMin = h([false;g]);
    leftMin(end)=[];
    rightMin(1)=[];
    
    % compute fwhh, max, and min  for every peak
    valMax = arrayfun(@(lm,rm) max(y(lm:rm)),leftMin,rightMin);
    posPeak = arrayfun(@(lm,rm,vm) find(vm==y(lm:rm),1)+lm-1,leftMin,rightMin,valMax);
    lfwhh = arrayfun(@(lm,vm,pp) interp1q(y(lm:pp),X(lm:pp),vm/2),leftMin,valMax,posPeak);
    rfwhh = arrayfun(@(rm,vm,pp) interp1q(y(rm:-1:pp),X(rm:-1:pp),vm/2),rightMin,valMax,posPeak);
    lfwhh(isnan(lfwhh)) = X(leftMin(isnan(lfwhh)));
    rfwhh(isnan(rfwhh)) = X(rightMin(isnan(rfwhh)));
    
    j = NaN;
    while ~isempty(j) %until no more oversegmented peaks
        % compute threshold for centroid for each peak
        peakThld = valMax * peakLocationRatio - sqrt(eps);
        % calculate the centroids for each peak
        pkX = arrayfun(@(lm,rm,th) (((y(lm:rm)).*(y(lm:rm)>=th))'*X(lm:rm))/sum((y(lm:rm)).*(y(lm:rm)>=th)),leftMin,rightMin,peakThld);
        % look for potential oversegmented peaks
        dpkX = [inf;diff(pkX);inf];
        j = find((dpkX(2:end-1)<=oversegFilter)  &...
                 ((dpkX(2:end-1)<=dpkX(1:end-2)))  & ...
                 ((dpkX(2:end-1)<dpkX(3:end)) ));
        leftMin(j+1) = []; 
        rightMin(j) = [];
        lfwhh(j+1) = [];
        rfwhh(j) = [];
        if numel(j)==1
            valMax(j) = max(valMax([j j+1]));
        else
            valMax(j) = max(valMax([j j+1]),[],2);
        end
        valMax(j+1) = [];
    end
    
    % remove peaks below the height and fwhh thresholds
    k = valMax>=heightFilter & ~((rfwhh-lfwhh)<fwhhFilter);
    % add peak list to output
    P{i} = [pkX(k) valMax(k)];
    if nargout>1
        PFWHH{i} = [lfwhh(k) rfwhh(k)];
    end
    if nargout>2
       PEXT{i} = [X(leftMin(k)) X(rightMin(k))];
    end
       
    if (i == plotId)
        figure
        hold on
        switch style
            case 'peak'
               plot(X,y_n,'displayname','Original signal');
               if waveletEnable
                  plot(X,y,'g','displayname','Denoised signal');
               end
               plot(pkX(k),valMax(k),'xr','linewidth',2,'displayname','Peaks');
            case 'exttriangle'
               xx = [X(leftMin(k)) pkX(k) X(rightMin(k)) nan(sum(k),1)]';
               yy = [y(leftMin(k)) valMax(k) y(rightMin(k)) nan(sum(k),1)]';
               plot(X,y_n,'displayname','Original signal');
               if waveletEnable
                  plot(X,y,'g','displayname','Denoised signal');
               end
               plot(xx(:),yy(:),'r','linewidth',2,'displayname','Peaks');
            case 'fwhhtriangle'
               xx = [lfwhh(k) pkX(k) rfwhh(k) nan(sum(k),1)]';
               yy = [max(interp1q(X,y,lfwhh(k)),valMax(k)/2) valMax(k) max(interp1q(X,y,rfwhh(k)),valMax(k)/2) nan(sum(k),1)]';
               plot(X,y_n,'displayname','Original signal');
               if waveletEnable
                  plot(X,y,'g','displayname','Denoised signal');
               end
               plot(xx(:),yy(:),'r','linewidth',2,'displayname','Peaks');
            case 'extline'
               xx = [X(leftMin(k)) X(leftMin(k)) nan(sum(k),1) X(rightMin(k)) X(rightMin(k)) nan(sum(k),1)]';
               yy = [zeros(sum(k),1) y(leftMin(k)) nan(sum(k),1) zeros(sum(k),1) y(rightMin(k)) nan(sum(k),1)]';
               hlin = plot(xx(:),yy(:),'r','linewidth',2,'displayname','Extents');
               setappdata(hlin,'legend_hgbehavior',0)
               plot(X,y_n,'displayname','Original signal');
               if waveletEnable
                  plot(X,y,'g','displayname','Denoised signal');
               end
               plot(pkX(k),valMax(k),'xr','linewidth',2,'displayname','Peaks'); 
            case 'fwhhline'
               xx = [lfwhh(k) lfwhh(k) rfwhh(k) rfwhh(k) nan(sum(k),1)]';
               yy = [interp1q(X,y,lfwhh(k)) valMax(k)/2 valMax(k)/2 interp1q(X,y,rfwhh(k)) nan(sum(k),1)]';
               hlin = plot(xx(:),yy(:),'r','linewidth',2,'displayname','FWHH');
               setappdata(hlin,'legend_hgbehavior',0)
               plot(X,y_n,'displayname','Original signal');
               if waveletEnable
                  plot(X,y,'g','displayname','Denoised signal');
               end
               plot(pkX(k),valMax(k),'xr','linewidth',2,'displayname','Peaks');                            
        end               
        title(sprintf('Signal ID: %d',i));
        xlabel('Separation Units')
        ylabel('Relative Intensity')
        legend('show')
        axis([min(X) max(X) min(y) max(y_n)])
        setAllowAxesRotate(rotate3d(gcf),gca,false)
        grid on
        hold off
    end
    
end % if nargout>0 || (i == plotId)    
end % for i = 1:numSignals 
    
if numSignals==1
    P = P{1};
    if nargout>1
        PFWHH = PFWHH{1};
    end
    if nargout>2
        PEXT = PEXT{1};
    end
end
   

if nargout == 0 
    clear P
end

