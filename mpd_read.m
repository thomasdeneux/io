function [data sfr] = mpd_read(fname,header)
% function [data sfr] = mpd_read(fname[,header])

if nargin==0, fname = fn_getfile('*.MPD'); end
if nargin<2, header = mpd_headers(fname); end

nheader = 20480;

% for the moment we suppose that
% - there is only one acquisition
% - there are 2 channels
% -> check the file length

nx = header.nx;
ny = header.ny;
nfr = header.nfr;
nelem = header.nelem;

% read data
fid = fopen(fname,'r');
fseek(fid,nheader,'bof');
x = fread(fid,[nx nelem/nx],'uint16');
fclose(fid);

% rearrange data: frames acquired as pairs!
data = zeros(nx,ny,nfr);
sfr  = zeros(nx,ny,nfr);
xpos = 0;
for i=0:floor(nfr/2)-1
    data(:,:,2*i+1) = x(:,xpos+0*ny+(1:ny));
    data(:,:,2*i+2) = x(:,xpos+1*ny+(1:ny));
    sfr(:,:,2*i+1)  = x(:,xpos+2*ny+(1:ny));
    sfr(:,:,2*i+2)  = x(:,xpos+3*ny+(1:ny));
    xpos = xpos + 4*ny;
    if mod(i,8)==7, xpos = xpos+8; end % aie aie aie
end
if mod(nfr,2)
    data(:,:,nfr) = x(:,:,xpos+0*ny,(1:ny));
    sfr(:,:,nfr)  = x(:,:,xpos+1*ny,(1:ny));
end
