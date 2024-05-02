function a = pv_read(f,channel)
% function a = pv_read(f[,channel])
%---
% reads 2-photon files from the XML format

% which file?
if nargin<1
    f = fn_getfile('*.xml');
end

folder = fileparts(f);
swd = pwd;
cd(folder)

% read XML headers
header = tps_headers(f);

% which channels?
if nargin>1
    channels = channel;
else
    channels = [];
    for j=1:size(header.files,2)
        if ~isempty(header.files{1,j}), channels = [channels j]; end
    end
end
nc = length(channels);

% read data
a = zeros(header.nx,header.ny,header.nfr,nc);

for i=1:header.nfr
    for j=1:nc
        a(:,:,i,j) = double(imread(header.files{i,channels(j)}));
    end
end

cd(swd)