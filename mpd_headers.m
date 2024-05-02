function header = mpd_headers(fname)
% function header = mpd_headers(fname)

if nargin==0, fname = fn_getfile('*.MPD'); end

nheader = 20480;
nhplus  = 8192; % at the end of each acquisition

% for the moment we suppose that
% - there is only one acquisition
% - there are 2 channels
% -> check the file length

d = dir(fname);
nelem = (d.bytes - (nheader+nhplus))/2; % uint16 is 2 bytes

if mod(nelem,512*512*2)==0
    % are images 512*512?
    nx = 512;
    ny = 512;
    nfr = nelem/(nx*ny*2);
elseif mod(nelem,256*256)==0
    % are images 256*256?
    nx = 256;
    ny = 256;
    nfr = nelem/(nx*ny*2);
else
    answer = inputdlg({'nx' 'ny' 'nfr'},'Enter image size',1,{'256' '256' '1'});
    nx = str2double(answer{1});
    ny = str2double(answer{2});
    nfr = str2double(answer{3});
end

header = struct('nx',nx,'ny',ny,'nfr',nfr,'nelem',nelem);

