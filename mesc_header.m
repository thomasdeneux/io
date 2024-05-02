function [h sessions] = mesc_header(file)
% function [h sessions] = mesc_header(file)
%---
% Returns meta-data extracted from .mesc file.
% If no output is requested, a small summary is printed only.
% 
% Input:
% - file        file name
% 
% Output:
% - h           full header
% - sessions    small summary


if nargin==0
    help mesc_header
    return
end

matfile = [file ' - headers.mat']; % output is saved into a mat file to save scanning the file multiple times
if exist(matfile,'file') 
    [h sessions] = fn_loadvar(matfile);
else
    
    tzoffset=1; % time shift with respect to GMT - TODO: use this offset to correct time!
    
    % open HDF5 file for reading
    try
        fileID=H5F.open(file, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    catch
        error('Unable to open file "%s".', file)
    end
    
    h = getobject(fileID,'/');
    nsession = double(h.VecMSessionsSize);
    sessions = cell(1,nsession);
    unitfields = {}; % all field names of the Units structure 
    for i=1:nsession
        hi = getobject(fileID,['MSession_' num2str(i-1)]);
        nunit = double(hi.VecMUnitsSize);
        sessi = cell(1,nunit);
        for j=1:nunit
            hij = getobject(fileID,['MSession_' num2str(i-1) '/MUnit_' num2str(j-1)]);
            try
                hi.Units(j) = hij;
            catch
                % dissimilar structures
                F0 = fieldnames(hi.Units);
                if ~isempty(F0)
                    % make sure that hi.Units(j) is created (would not be
                    % the case if hij has no field) 
                    hi.Units(j).(F0{1}) = [];
                end
                F = fieldnames(hij);
                for k=1:length(F)
                    f = F{k};
                    hi.Units(j).(f) = hij.(f);
                end
            end
            try
                sessi{j} = sprintf('Session%i/Unit%i %ix%ix%i %s',i-1,j-1,hij.XDim,hij.YDim,hij.ZDim,char(hij.Comment'));
            catch
                sessi{j} = sprintf('Session%i/Unit%i ERROR',i-1,j-1);
            end
        end
        h.sessions(i) = hi;
        sessions{i} = char(sessi);
        if nunit>0, unitfields = union(unitfields,fieldnames(hi.Units),'stable'); end
    end
    % add fields if necessary to make all 'Units' structure similar
    for i=1:nsession
        hi = h.sessions(i); 
        if isempty(hi.Units), continue, end
        F = fieldnames(hi.Units);
        if ~isequal(F,unitfields)
            F1 = setdiff(unitfields,F);
            for k=1:length(F1)
                hi.Units(1).(F1{k}) = [];
            end
            hi.Units = orderfields(hi.Units,unitfields);
        end
        h.sessions(i) = hi;
    end
    sessions = char(sessions);
    
    % save
    fn_savevar(matfile,h,sessions);
end

if nargout==0
    clear h
    disp(sessions)
end

%---
function obj = getobject(fileID,objname)


% open HDF5 group or dataset
try
    objID=H5G.open(fileID, objname);
catch
    try
        objID=H5D.open(fileID, objname);
    catch
        fprintf('Unable to open HDF5 object "%s".\n', objname)
        obj = struct;
        return
    end
end

% iterate over all contained attributes
obj = struct;
numAttrs=H5A.get_num_attrs(objID);
for attrIdx=0:numAttrs-1
    try
        % read data from the currently inspected attribute
        attribID=H5A.open_idx(objID, attrIdx);
        attribName=H5A.get_name(attribID);
        attribValue=H5A.read(attribID, 'H5ML_DEFAULT');
        H5A.close(attribID);

        % print attribute name, class, and value
        obj.(attribName) = attribValue;
    catch
        % issue warning for unreadable attribute (due to some HDF5
        % forward incompatibility)
        
        %fprintf('attribute %d of object %s unreadable\n', attrIdx, objname)
    end
end

% HDF5 objects and files are closed automatically on function exit
