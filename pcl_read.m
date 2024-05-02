function [movie, param, time, position] = pcl_read(fname, param)
% function [movie, param, time, position] = pcl_read(fname, param)
% function [movie, param, time, position] = pcl_read(fname, 'all')

if nargin==0
    fname = fn_getfile('*.pcl');
end

%% Read file and estimate finest time resolution

persistent last_file last_data

fmat = strrep(fname,'.pcl','.mat');
if strcmp(fmat, last_file)
    disp 'raw data already in memory'
    [time, position, dt0] = deal(last_data{:});
else
    if exist(fmat,'file')
        disp 'reading matlab file'
        [time, position, dt0] = fn_loadvar(fmat,'time','position','dt0');
        if isa(time, 'double')
            disp 're-save raw data to Matlab file of smaller size'
            time = single(time);
            position = uint16(position);
            save(fmat,'time','position','dt0','-v7.3')
        end
    else
        %%
        d = dir(fname);
        nrec = d.bytes / 14;

        %%
        fid = fopen(fname,'r');
        %%
        time = zeros(nrec,1);
        x = zeros(nrec,1,'uint8');
        y = zeros(nrec,1,'uint8');
        %%
        disp 'reading times'
        fseek(fid,0,'bof');
        time = fread(fid,Inf,'double=>single',6,'b');
        disp 'reading positions'
        fseek(fid,8,'bof');
        position = fread(fid,[2 Inf],'2*uint16=>uint16',10,'b');
        %%
        time = time - time(1);
        %%
        disp 'determine finest temporal resolution'
        dd = diff(time); dd(dd==0 | dd>.1) = [];
        dt0 = median(dd);
        dt0 = dichotomy(@(dt)check_time_resolution(time,dt),dt0*.9,dt0*1.1,1e-7);
        fprintf('\b -> %fms\n',dt0*1e3)
        %% (save)
        disp 'save raw data to Matlab file'
        save(fmat,'time','position','dt0','-v7.3')
    end
        
    last_file = fmat;
    last_data = {time position dt0};
end

answer = questdlg('Keep raw data in memory for a later usage?', '', 'Yes', 'No', 'Yes');
if ~strcmp(answer, 'Yes')
    last_file = [];
    last_data = [];
end

%%
% ask user for certain parameters
if nargin<2   
    param = fn_structedit( ...
        'dt',   {dt0     'double'    'Temporal precision (seconds)'}, ...
        'tstart',   {0  'double'    'Start time (minutes)'}, ...
        'tstop',    {ceil(max(time)/60) 'double'  'Stop time (minutes)'}, ...
        'nx',   {max(position(1,:))+1 'double'    'image width'}, ...
        'ny',   {max(position(2,:))+1 'double'    'image height'}, ...
        'xbin', {3      'double'    'spatial binning'}, ...
        'dosparse', {false 'logical' 'Use sparse array'});
end
nx0 = double(param.nx);
ny0 = double(param.ny);
dt = param.dt;
tstart = param.tstart*60;
tstop = param.tstop * 60;
nt = ceil((tstop-tstart)/dt);
xbin = param.xbin;
nx = 1 + ceil(nx0/xbin);
ny = 1 + ceil(ny0/xbin);
%%
% remove times out of the window
ok = (time>=tstart & time<tstop & position(2,:)'>1);
time = time(ok);
time = time - time(1);
position = position(:,ok);
subs = {1+floor(position(1,:)/xbin)' 1+floor(position(2,:)/xbin)' ...
    1+floor(time/dt)};
clear position
if param.dosparse
    disp 'making ndSparse movie'
    movie = ndSparse.accumarray(subs, 1, [nx ny nt]);
else
    disp 'making uint8 movie'
    movie = accumarray(subs, uint8(1), [nx ny nt], @(x) sum(x,'native'));
    % check 
    m = max(movie(:));
    if m == uint8(Inf)
        % we need to use uint16!
        disp 'making uint16 movie'
        movie = accumarray(subs, uint16(1), [nx ny nt], @(x) sum(x,'native'));
        m = max(movie(:));
        if m == uint16(Inf)
            % we need to use uint32!
            disp 'making uint32 movie'
            movie = accumarray(subs, uint32(1), [nx ny nt], @(x) sum(x,'native'));
        end
    end
end

%---
function e = check_time_resolution(time,dt)


t = mod(time/dt+.5,1);
d = abs(diff(t));
d(d>.5) = [];
e = mean(d);
% figure(1), plot(t)


