function s = oi_listmat(fname)
% function s = oi_listmat([fname])
% function s = oi_listmat(dirname)
%---
% Returns a structure with information on a list of block (mat) files.
% 
% Input:
% - fname   list of file names (can be either 2D char array or cell array
%           of strings); prompt user for files if not specified
% - dirname folder name: scans all files inside it
%
% Output:
% - s       structure with fields:
%   .datadir        path to data
%   .files          file names, sorted according to acquisition time
%   .nx, .ny, .nt   dimension
%   .ncond          number of conditions
%   .exp            structure with fields 'k' (experiment number) and
%                   'idx' (indices of file names concerned with experiment
%                   k) 
%   .cond           condition number for each file (case where each file
%                   contains only one trial)
%
%---
% Thomas Deneux
% note that this function is still under development and is not fully
% operational yet

s = struct;

% get list of files (both as cell array and char array), detect data directory
if nargin==0
    fname2 = fn_getfile('*.mat');
    if ~fname2, return, end
    s.datadir = fileparts(fname2(1,:));
    fname2(:,1:length(s.datadir)+1) = []; % remove path component
    fname = cellstr(fname2);
elseif ~iscell(fname) && isvector(fname) && isdir(fname)
    swd = pwd; cd(fname), s.datadir = pwd; cd(swd)
    fname = dir([s.datadir '/*.mat']);
    if isempty(fname), disp('could not find any file'), return, end
    fname = {fname.name};
    fname2 = strvcat(fname{:});  
else
    if iscell(fname), fname2 = strvcat(fname{:}); end 
    p = fileparts(fname2(1,:));
    if isempty(p)
        s.datadir = fileparts(which(deblank(fname(1,:))));
    else
        swd = pwd; cd(p), s.datadir = pwd; cd(swd)
        fname2(:,1:length(p)+1) = [];
    end
    fname = cellstr(fname2);
end

% scan the first name for position of flags
f = fname{1};
spos = 0;   % successflag
if ismember(f(1),'TF'), spos = 1; end
flags = 'CEB';
xpos = zeros(4,2);
tpos = [0 0];
for i=1:3
    [idx tmp] = regexp(f,[flags(i) '(\d+)'],'start','tokens');
    if isempty(tmp), continue, end
    if ~isscalar(tmp) || ~isscalar(tmp{1}), error('cannot read file name'), end
    len = length(tmp{1}{1});
    while len>4 % shit, includes date and/or time!!
        len = len-6; % not sure that this works...
        if len<2, error('cannot read file name'), end
        tpos = [idx+len tpos(2)+6];
    end
    xpos(i,:) = [idx len];
end

% remove bad trials
if spos
    badidx = (fname2(:,1)=='F');
    fname2(badidx,:) = [];
    fname(badidx) = [];
end


% exp & blocks
exp = str2num(fname2(:,xpos(2,1)+(1:xpos(2,2))))'; %#ok<*ST2NM>
blocks = str2num(fname2(:,xpos(3,1)+(1:xpos(3,2))))';
[tmp ord] = sortrows([exp' blocks']);
exp = tmp(:,1)';
fname  = fname(ord);
fname2 = fname2(ord,:);
s.files = fname2;
s.nf = size(fname2,1);
exps = unique(exp);
if isscalar(exps)
    s.exp = exps;
else
    s.exp = struct('k',num2cell(exps));
    for i=1:length(exps)
        s.exp(i).idx = find(exp==exps(i));
    end
end

% cond
if xpos(2,1)
    cond = str2num(fname2(:,xpos(1,1)+(1:xpos(1,2))))'; 
    conds = unique(cond);
    s.ncond = length(conds);
    s.condvdaq = cond;
    cond2idx = zeros(1,max(conds+1));
    cond2idx(conds+1) = 1:s.ncond;
    s.cond = cond2idx(cond+1);
end
    
% time
if tpos(1)
    if tpos(2)==12
        tformat = 'ddmmyyHHMMSS';
    else
        disp('please help me to decide date/time format')
        keyboard
    end
    tt = datenum(s.files(:,tpos(1)+(1:tpos(2))),tformat);
    if any(diff(floor(tt))), disp('several dates! please help me'), keyboard, end
    s.datestr = datestr(floor(tt(1)));
    s.timestr = cellstr(datestr(mod(tt,1)));
    s.time = tt';
end










