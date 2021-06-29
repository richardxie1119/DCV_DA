function [manualLbl] = drawAndAvg(figNumIn, xyData, mzs, intens, figNumOut, manualLbl)
%figNumIn: figure to draw on
%xyData: data plotted in figNumIn
%mzs: x axis of average plot
%intens: intensity of average plot
%figNumOut: where to plot average, 0 to make a new plot
    
    if isempty(manualLbl)
       manualLbl =  zeros(size(xyData,1),1);
    end
    curInd = max(manualLbl)+1;
    
    figure(figNumIn);
    roi = impoly;
    tempPoly = roi.getPosition();
    inpoly = inpolygon(xyData(:,1), xyData(:,2), tempPoly(:,1), tempPoly(:,2));
    
    manualLbl(inpoly) = curInd;

    %show average
    if figNumOut
        figure(figNumOut);
        hold on;
    else
        figure; hold on;
    end
    if sum(inpoly) > 1
        m = reshape(repmat(mzs,3,1), 1,[])+0.01*curInd;
        i = zeros(1,length(mzs)*3);
        i(2:3:end) = mean(intens(:,inpoly),2);
        plot(m,i);
        xy = mean(xyData(inpoly,:));
        figure(figNumIn); hold on;
        text(xy(1), xy(2), num2str(curInd));
    else
        m = reshape(repmat(mzs,3,1), 1,[])+0.01*curInd;
        i = zeros(1,length(mzs)*3);
        i(2:3:end) = intens(:,inpoly);
        plot(m,i);
        xy = xyData(inpoly,:);
        figure(figNumIn); hold on;
        text(xy(1), xy(2), num2str(curInd));
    end

end