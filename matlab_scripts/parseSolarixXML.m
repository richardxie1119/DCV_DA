function [ mz, intens, res ] = parseSolarixXML( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %initialize return values
    mz = [];
    intens = [];
    res = [];
    %open file
    fid = fopen(filename, 'r');
    %for each line
    line = fgetl(fid);
    while(~isempty(line))
        %check if line has a peak
        if(line(1:3) == '<pk')
            toks = strsplit(line, ' ');
            %parse out the needed information
            mz = [mz str2double(toks{2}(5:end-1))];
            intens = [intens str2double(toks{3}(4:end-1))];
            res = [res str2double(toks{7}(6:end-1))];
        elseif (strcmp(line, '</ms_peaks>'))
            break
        end      
        
        line = fgetl(fid);
    end
    fclose(fid);
end

