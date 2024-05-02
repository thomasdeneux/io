function data = mesc_read(file, msessionIdx, munitIdx, channelIdx, frameIdx)
% function data = mesc_read(file, msessionIdx, munitIdx, channelIdx[, frameIdx])
%---
% returns frames from a movie in a .mesc file.
%
%   Parameters:
%       file        .mesc file
%       MSESSIONIDX measurement session index (usually 0)
%       MUNITIDX    measurement unit index (indexed from 0)
%       CHANNELIDX  channel index (indexed from 0)
%       FRAMEIDX    index of single or consecutive frame(s) to return (indexed from 0, [] for all frames)

if nargin==0, help mesc_read; end
if nargin<5, frameIdx=[]; end

% open HDF5 file for reading
try
    fileID=H5F.open(file, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
catch
    error('Unable to open file "%s".', file)
end

% open the HDF5 group representing the measurement unit
groupname=sprintf('/MSession_%d/MUnit_%d', msessionIdx, munitIdx);
try
    groupID=H5G.open(fileID, groupname);
catch
    error('Unable to open measurement unit.')
end

% open the HDF5 dataset holding the channel contents
try
    datasetID=H5D.open(fileID, ...
        sprintf('%s/Channel_%d', groupname, channelIdx));
catch
    error('Unable to open channel.')
end

% get the necessary attributes
attribID=H5A.open_name(groupID, 'XDim');
xdim=H5A.read(attribID, 'H5ML_DEFAULT');
H5A.close(attribID);
attribID=H5A.open_name(groupID, 'YDim');
ydim=H5A.read(attribID, 'H5ML_DEFAULT');
H5A.close(attribID);
attribID=H5A.open_name(groupID, 'ZDim');
zdim=H5A.read(attribID, 'H5ML_DEFAULT');
H5A.close(attribID);

if isempty(frameIdx)
    framestart = 0;
    framecount = zdim;
elseif ~all(diff(frameIdx)==1)
    error 'frame indices must be consecutives'
else
    framestart = frameIdx(1);
    framecount = length(frameIdx);
    if framestart+framecount > zdim
        error('Section index is too large; section does not exist.')
    end
end

stride = [1 1 1];
count = [1 1 1];

% read blocks of 1GB
nbyteperframe = ydim*xdim*2;
blocksize = 2^30;
nframeperblock = floor(blocksize/double(nbyteperframe));
nblock = ceil(double(framecount)/nframeperblock);
if nblock>1
    % pre-allocate
    data = zeros(xdim,ydim,framecount,'uint16'); 
    fn_progress('reading data per blocks',nblock)
end
for iblock = 1:nblock
    if nblock>1, fn_progress(iblock), end
    blockstart = double([framestart+(iblock-1)*nframeperblock 0 0]);
    blockcount = fn_switch(iblock<nblock,nframeperblock,framecount-(nblock-1)*nframeperblock);
    try
        block = double([blockcount ydim xdim]);
        dataspaceID = H5D.get_space(datasetID);
        H5S.select_hyperslab(dataspaceID,'H5S_SELECT_SET',blockstart,stride,count, block);
        dataspaceID_memory = H5S.create_simple(3,block,[]);
        datak = H5D.read(datasetID, 'H5ML_DEFAULT', dataspaceID_memory, ...
            dataspaceID, 'H5P_DEFAULT');
        datak = 65535 - datak;
        datak = datak(:,end:-1:1,:);
        if nblock==1, data = datak; else data(:,:,blockstart(1)+(1:blockcount)) = datak; end
    catch ME,
        disp 'reading full data block failed, error was:'
        disp(ME.message)
        disp 'try reading frame by frame'
        
        if nblock==1, data = zeros(xdim,ydim,blockcount,'uint16'); end % data was not allocated yet
        for k=1:blockcount
            startk = double([blockstart(1)+k-1 0 0]);
            stride = [1 1 1];
            count = [1 1 1];
            blockk = double([1 ydim xdim]);
            dataspaceID = H5D.get_space(datasetID);
            H5S.select_hyperslab(dataspaceID,'H5S_SELECT_SET',startk,stride,count, blockk);
            dataspaceID_memory = H5S.create_simple(3,blockk,[]);
            try
                datak = H5D.read(datasetID, 'H5ML_DEFAULT', dataspaceID_memory, ...
                    dataspaceID, 'H5P_DEFAULT');
                datak = 65535 - datak(:,end:-1:1);
                data(:,:,blockstart(1)+k) = datak;
            catch
                disp(['failure for frame ' num2str(blockstart(1)+k) ', replacing by previous'])
                if startk(1)==0, error 'ca suffit', end
                data(:,:,startk(1)+1) = data(:,:,startk(1));
            end
        end
    end
end


% HDF5 objects and files are closed automatically on function exit
