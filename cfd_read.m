function Y = cfd_read(fname)
% function Y = cfd_read(fname)

% Input
if nargin<1
    fname = fn_getfile('*.CFD');
end

% Multiple files
if size(fname,1)>1
    nf = size(fname,1);
    a = cfd_read(fname(1,:));
    [nx ny nt nc] = size(a);
    Y = zeros(nx,ny,nt,nc,nf);
    Y(:,:,:,:,1) = a;
    for i=2:nf
        a = cfd_read(fname(i,:));
        Y(:,:,:,:,i) = a;
    end
    return
end

% Read headers
headers = cfd_headers(deblank(fname));

% Open file
fid = fopen(deblank(fname),'r');
fseek(fid,headers.CFH.hdrsize,-1);

% Read data
nfrm = headers.ACV.num.images;
% disp('warning: don''t know where to find the number of channels in headers - assume there are two')
% VERY STRANGE: 1 means 1 channel, and 3 means 2 channels!?
nchans = (headers.ACV.num.chans+1)/2; 
nx = headers.ACV.pixels.x;
ny = headers.ACV.pixels.y;
Y = fread(fid,nx*ny*nchans*nfrm,'uint8=>uint8');
if nfrm~=length(Y)/(nx*ny*nchans)
    disp(['warning: file truncated: ' fname])
    nfrm = length(Y)/(nx*ny*nchans);
end
Y = reshape(Y,[nchans nx ny nfrm]);
Y = permute(Y,[2 3 4 1]);
%Y = single(Y);

% Close file
fclose(fid);