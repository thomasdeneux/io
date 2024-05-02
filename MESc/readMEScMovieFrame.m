function data = readMEScMovieFrame(...
    path, msessionIdx, munitIdx, channelIdx, frameIdx)
% READMESCMOVIEFRAME returns a frame from a movie in a .mesc file.
%   readMEScMovieFrame( path, msessionIdx, munitIdx, channelIdx, frameIdx )
%   returns a frame from a .mesc movie measurement unit.
%   Parameters:
%       PATH        .mesc file path
%       MSESSIONIDX measurement session index (usually 0)
%       MUNITIDX    measurement unit index (indexed from 0)
%       CHANNELIDX  channel index (indexed from 0)
%       FRAMEIDX    index of the frame to return (indexed from 0)

if nargin < 5
    error('usage: data = readMEScMovieFrame(path, msessionIdx, munitIdx, channelIdx, frameIdx)')
end

% open HDF5 file for reading
try
    fileID=H5F.open(path, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
catch
    error('Unable to open file "%s".', path)
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
if frameIdx >= zdim
    error('Section index is too large; section does not exist.')
end

start = [frameIdx 0 0];
stride = [1 1 1];
count = [1 1 1];
block = [1 ydim xdim];
dataspaceID = H5D.get_space(datasetID);
H5S.select_hyperslab(dataspaceID,'H5S_SELECT_SET',start,stride,count, block);
dataspaceID_memory = H5S.create_simple(3,block,[]);

data = 65535 - permute( ...
    H5D.read(datasetID, 'H5ML_DEFAULT', dataspaceID_memory, ...
    dataspaceID, 'H5P_DEFAULT'),[2 1]);

% HDF5 objects and files are closed automatically on function exit
