function [ colormap ] = myhsv( numvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    maps = {[68,119,170]/255
        
        [68,119,170
        204,102,119]/255
        
        [68,119,170
        221,204,119
        204,102,119]/255
        
        [68,119,170
        17,119,51
        221,204,119
        204,102,119]/255
        
        [51,34,136
        136,204,238
        17,119,51
        221,204,119
        204,102,119]/255
        
        [51,34,136
        136,204,238
        17,119,51
        221,204,119
        204,102,119
        170,68,153]/255
        
        [51,34,136
        136,204,238
        68,170,153
        17,119,51
        221,204,119
        204,102,119
        170,68,153]/255
        
        [51,34,136
        136,204,238
        68,170,153
        17,119,51
        153,153,51
        221,204,119
        204,102,119
        170,68,153]/255
        
        [51,34,136
        136,204,238
        68,170,153
        17,119,51
        153,153,51
        221,204,119
        204,102,119
        136,34,85
        170,68,153]/255
        
        [51,34,136
        136,204,238
        68,170,153
        17,119,51
        153,153,51
        221,204,119
        102,17,0
        204,102,119
        136,34,85
        170,68,153]/255
        
        [51,34,136
        102,153,204
        136,204,238
        68,170,153
        17,119,51
        153,153,51
        221,204,119
        102,17,0
        204,102,119
        136,34,85
        170,68,153]/255
        
        [51,34,136
        102,153,204
        136,204,238
        68,170,153
        17,119,51
        153,153,51
        221,204,119
        102,17,0
        204,102,119
        170,68,102
        136,34,85
        170,68,153]/255
        
        };
    colormap = maps{numvals};
end
