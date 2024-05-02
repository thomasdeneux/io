function T = oi_binblocks(arg)
% function T = oi_binblocks(files)
% function T = oi_binblocks(T)
%---
% opens a dialog: select block file to bin, then select a saving name, that
% will be used for both a new .tptrial file and a new folder containing the
% binned data

% Thomas Deneux
% Copyright 2012-2012

dooutput = (nargout>0);
dosavetptrial = ~dooutput;

% prompt user for files and where to save the converted files
okT = nargin>0 && isa(arg,'tps_trial');
if okT
    T = arg;
    files = {T.file};
else 
    if nargin<1
        files = {};
        f = fn_getfile('*.BLK*','Select files to bin');
        if isequal(f,0), return, end
        while ~isequal(f,0)
            files = [files; cellstr(f)]; %#ok<AGROW>
            f = fn_getfile('*.BLK*','Select more files to bin');
        end
    else
        files = cellstr(arg);
    end
    if isscalar(files) && strcmp(fn_fileparts(files{1}, 'ext'),'tptrial')
        T = fn_loadvar(files{1});
        okT = true;
        files = {T.files};
    end
end
nf = length(files);

% binning value
s = fn_structedit('xbin',{2 'stepper 1 1 100' 'spatial binning'},'tbin',{1 'stepper 1 1 100' 'temporal binning'});
xbin = s.xbin;
tbin = s.tbin;
pause(.1) % what a bug!

% location for saving
if dosavetptrial
    msg = 'Select name for folder and .tptrial file';
    fsave = fn_savefile('*.tptrial',msg);
else
    msg = 'Select folder where to save binned data';
    fsave = fn_getdir(msg);
end
if isequal(fsave,0), return, end
[p base] = fileparts(fsave);
fsave = [p '/' base];
if ~exist(fsave,'dir'), mkdir(fsave), end

% bin data and store into a tps_trial object, so that header information
% will be saved as well
if ~okT, T = tps_trial.empty(1,0); end
fn_progress('binning file',nf)
for k=1:nf
    fn_progress(k)
    [dum base ext] = fileparts(files{k}); %#ok<ASGLU>
    fsavek = [fsave '/' base ext '.bin' num2str(xbin) '.mat'];
    if okT, Tk = T(k); else Tk = tps_trial(files{k}); end
    if exist(fsavek,'file')
        answer = fn_dialog_questandmem('binned file exist, use existing file?');
        if ~answer, return, end
        v = load(fsavek);
        a = v.data;
    else
        a = Tk.data;
        a = fn_float(a);
        a = fn_bin(a,[xbin xbin tbin 1]);
    end
    T(k) = tps_trial(a,Tk);
    T(k).dx = Tk.dx*xbin;
    T(k).dy = Tk.dy*xbin;
    T(k).dt = Tk.dt*tbin;
    clear Tk
    fn_progress('saving')
    if exist(fsavek,'file'), T(k).file = fsavek; else savedata(T(k),fsavek), end
end

% merge the stim tables
mergestimtables(T)

% save
if dosavetptrial
    save([fsave '.tptrial'],'T','-MAT')
end

% output?
if ~dooutput, clear T, end