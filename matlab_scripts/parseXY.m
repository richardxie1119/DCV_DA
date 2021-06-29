function [ x, y] = parseXY( instr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    toks = strsplit(instr, '_');
    x = str2double(toks{2}(1:end-1));
    y = str2double(toks{3});
end

