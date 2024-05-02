function [h sessions] = mes_header(fmes)
% function [h sessions] = mes_header(file)
%---
% Returns meta-data extracted from .mes file.
% If no output is requested, a small summary is printed only.
% 
% Input:
% - file        file name
% 
% Output:
% - h           full header
% - sessions    small summary


if nargin==0
    help mes_header
    return
end

matfile = [fmes ' - headers.mat']; % output is saved into a mat file to save scanning the file multiple times
if exist(matfile,'file') 
    [h sessions] = fn_loadvar(matfile);
else
    % load header info
    measurements = load(fmes,'D*','-MAT');
    measnames = fieldnames(measurements);
    h = struct2cell(measurements);
    nmes = length(h);
    % build short descriptions (measurement number, size, comment)
    sessions = cell(1,nmes);
    for k=1:nmes
        xk = h{k};
        xk1 = xk(1);
        siz = [xk1.Width xk1.Height];
        if length(xk)>1, siz = [siz length(xk)]; end %#ok<AGROW>
        siz = fn_strcat(siz,'x');
        comment = fn_switch(isempty(xk1.Comment),'',[' - ' xk1.Comment]);
        sessions{k} = [measnames{k} ' - ' siz ' ' xk1.Context comment];
    end    
    % save
    fn_savevar(matfile,h,sessions)
end
            