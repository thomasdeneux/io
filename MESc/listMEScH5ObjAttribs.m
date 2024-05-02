function listMEScH5ObjAttribs(path, objname, tzoffset)
% LISTMESCH5OBJATTRIBS lists the attributes of a HDF5 object in a MESc file.
%   listMEScH5ObjAttribs( path, objname, tzoffset ) lists
%   the names, classes, and values of each HDF5 attribute
%   in the HDF5 group or dataset OBJNAME in the HDF5 file at PATH.
%   Optional argument TZOFFSET denotes the hour offset of the local
%   timezone relative to UTC (e.g., 1 for Central European Time and 2 for
%   Central European Summer Time). The default value of TZOFFSET is 0.

if nargin < 2
    error('usage: listMEScH5ObjAttribs(path, objname[, tzoffset])')
end

if ~exist('tzoffset', 'var')
    tzoffset=0; % set default tzoffset
end

% open HDF5 file for reading
try
    fileID=H5F.open(path, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
catch
    error('Unable to open file "%s".', path)
end

% open HDF5 group or dataset
try
    objID=H5G.open(fileID, objname);
catch
    try
        objID=H5D.open(fileID, objname);
    catch
        error('Unable to open HDF5 object "%s".', objname)
    end
end

% iterate over all contained attributes
numAttrs=H5A.get_num_attrs(objID);
for attrIdx=0:numAttrs-1
    try
        % read data from the currently inspected attribute
        attribID=H5A.open_idx(objID, attrIdx);
        attribName=H5A.get_name(attribID);
        attribValue=H5A.read(attribID, 'H5ML_DEFAULT');
        attribClass=class(attribValue);
        H5A.close(attribID);

        % print attribute name, class, and value
        disp(['--> attribute ' num2str(attrIdx) ': name: ''' attribName ...
            ''', class: ' attribClass ', value:' ])
        if ischar(attribValue)
            % pretty-print 8-bit string
            fprintf('%s\n\n',attribValue')
        elseif strcmp(attribClass, 'int16')
            % pretty-print 16-bit string
            fprintf('%s\n\n',char(attribValue'))
        elseif strcmp(attribClass, 'uint64') && ...
                ( ~isempty(regexp(attribName, 'Time$', 'once')) ...
                || ~isempty(regexp(attribName, 'DatePosix$', 'once')))
            % pretty-print date to 1 second precision
            fprintf('%s\n\n',datestr(datenum(1970, 1, 1, 0, 0, ...
                double(attribValue))+tzoffset/24))
        else
            % print any other type of data
            disp(attribValue')
        end
    catch
        % issue warning for unreadable attribute (due to some HDF5
        % forward incompatibility)
        fprintf('--> attribute %d: ***unreadable***\n\n', attrIdx)
    end
end

% HDF5 objects and files are closed automatically on function exit
