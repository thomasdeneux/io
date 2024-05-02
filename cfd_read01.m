function [Yavgtrials Yalltrials] = cfd_read01(fname)
% function [Yavgtrials Yalltrials] = cfd_read01(fname)
%---
% reads an even number of CFD files corresponding to alternating blank and
% stimulated conditions
% - 
% IMPLEMENTATION NOT FINISHED

% Input
if nargin<1
    fname = fn_getfile('*.CFD');
end

nf = size(fname,1);
if mod(nf,2)~=0, error('there should be an even number of files'), end
nf = nf/2;

% Read
fn_progress('trial',nf*2);
a = cfd_read(fname(1,:));
s = size(a); clear a
Yavgtrials = zeros([s(1:4) 2]);
if nargout>1, Yalltrials = zeros([s(1:4) 2*nf]); end
for i=1:nf
    fn_progress(2*i-1)
    for k=0:1
        a = cfd_read(fname(2*i-1+k,:));
        Yavgtrials(:,:,:,:,1+k) = Yavgtrials(:,:,:,:,1+k) + a/nf;
        if nargout>1, Yalltrials(:,:,:,:,2*i-1+k) = a; end
    end
end
        