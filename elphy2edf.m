function elphy2edf(elphyfile, edffile)
% function elphy2edf([elphyfile[, edffile]])
%---
% Convert an Elphy .DAT file to European Data Format .edf file 

% Thomas Deneux
% Copyright 2018

if nargin<1
    elphyfile = fn_getfile('*.DAT', 'Select Elphy .DAT file to convert');
end
if nargin<2
    edffile = strrep(elphyfile,'.DAT','.edf');
end


