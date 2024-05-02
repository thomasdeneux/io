function [header data channels darkframe] = neuroplex_read(fname,channelflag)
% function [header data channels darkframe] = neuroplex_read(fname)
% function channels = neuroplex_read(fname,'channels')

% see e.g. http://www.redshirtimaging.com/support/dfo.html

% Input
if nargin==0
    fname = fn_getfile('*.da');
end
if nargin>=2
    if ~strcmp(channelflag,'channels'), error argument, end
    channelonly = true;
else
    channelonly = false;
end
dodata = ~channelonly && (nargout>=2);
dochannel = channelonly || (nargout>=3);

% open file
fid = fopen(fname,'r');

% read header
head = fread(fid,2560,'int16');
header = struct;

header.NFrames = head(5);

% for Photodiode array (NeuroPDA and NeuroPDA-II)  
% header.NPixels = head(97);
% header.FrameIntervalMS = header.NPixels * head(4) / 20000; 

% for CCD/CMOS cameras
header.NColumns = head(385);
header.NRows = head(386);
header.FrameIntervalMS = head(389)/1000;
DividingFactor = head(391);
if header.FrameIntervalMS>=10
    header.FrameIntervalMS = header.FrameIntervalMS*DividingFactor;
elseif DividingFactor~=1 && fn_dodebug
    disp 'strange, i expected DividingFactor to be 1 if it is not meant to be used'
    keyboard
end
header.AcquisitionRatio = head(392);

% data
[nx ny nfr] = deal(header.NColumns,header.NRows,header.NFrames);
if dodata
    data = reshape(fread(fid,nfr*nx*ny,'*int16'),[nfr nx ny]);
    data = permute(data,[2 3 1]); % x*y*time
else
    fseek(fid,2*nfr*nx*ny,'cof');
end

% channels
if dochannel
    channels = fread(fid,[header.AcquisitionRatio*nfr 8],'*int16');
else
    fseek(fid,2*header.AcquisitionRatio*nfr,'cof');
end

% dark frame
if dodata
    darkframe = fread(fid,[nx ny],'*int16');
end

% % check that we reached end of file
% disp(ftell(fid))
% fseek(fid,0,'eof');
% disp(ftell(fid))

% close file
fclose(fid);

% output
if channelonly
    header = channels;
end
