function data = oi_loadBLK(varargin)
% function data = oi_loadBLK([fnames][,flags,...])
%---
% load one or several block files (.BLK) through calling get_head and
% loadsum
% 
% Inputs:
% - fnames      file names - prompted if not specified
% - flag        'av' to specify to compute only the mean
%               a number (n) for binning n blocks together [not implemented yet] 
%               'array' to put all conditions and trials in the same array
%               (otherwise output is a cell array)
%               'uint8','uint16','uint32','single','double' to set the
%               output data type
%               'condn' to load only condition n
%
% Outputs:
% - data        cell array: 2D cell (trial x condition) of 3D arrays (horizontal x
%               vertical x time)
%               or 5D array (if 'array' flag): hor x vert x time x cond x
%               trial


% Input
fname = []; flag = ''; doarray = false; type_out = ''; conds = [];
for k=1:nargin
    a = varargin{k};
    if ~ischar(a), error('arguments must be char'), end
    if ischar(a) && isempty(fname) 
        fname = a;
    else
        switch a
            case 'array'
                doarray = true;
            case {'uint8' 'uint16' 'uint32' 'single' 'double'}
                type_out = a;
            case {'av' 'avg' 'average'}
                flag = 'av';
            otherwise
                if findstr(a,'cond')
                    conds = str2double(a(5:end));
                else
                    error('unknown flag ''%s''',a)
                end
        end
    end
end
if isempty(fname)
    try
        fname = fn_getfile('*.BLK','Select block file');
    catch
        [fname pathname] = uigetfile('*.BLK','Select block file','MultiSelect','on');
        fname = fullfile(pathname,fname);
        fname = char(fname);
    end
end
if ~fname, data={}; return, end

nblocks = size(fname,1);

% READING HEADER
header = oi_headBLK(deblank(fname(1,:)));
[nstim ni nj nfr] = deal(header.nstimuli,header.framewidth,header.frameheight,header.nframesperstim);
%     oi_get_head(deblank(fname(1,:)),'nstimuli,framewidth,frameheight,nframesperstim');
[lenh framesize filesize] = deal(header.lenheader,header.framesize,header.filesize);
[dtype gain dc scalefactor] = deal(header.datatype,header.meanampgain,header.meanampdc,header.scalefactor);
%     oi_get_head(deblank(fname(1,:)),'lenheader,framesize,filesize,datatype,meanampgain,meanampdc,scalefactor');

% GET_HEAD
% DAT_UCHAR     (11)
% DAT_USHORT    (12)
% DAT_LONG      (13)
% DAT_FLOAT     (14)

% LOADSUM
%	DATA_TYPE: 0 for ushort  (16bit) (CCD, raw VDAQ, ORA) -- default
%		   1 for ulong   (32bit) (sum16)
%		   2 for float   (32bit) (sumsum, convert)
%		   3 for double  (64bit) (matlab)
%		   4 for signed char     (Fuji)
%		   5 for unsigned char   (---)
%		   6 for signed short	 (DYEDAQ)

switch dtype
    case 11
        nbytes = 1;
        datatype = 'uchar';
        if isempty(type_out), type_out = 'uint8'; end
    case 12
        nbytes = 2;
        datatype = 'ushort';
        if isempty(type_out), type_out = 'uint16'; end
    case 13
        %disp('warning: reading longs as ulongs')
        nbytes = 4;
        datatype = 'ulong';
        if isempty(type_out), type_out = 'uint32'; end
    case 14
        nbytes = 4;
        datatype = 'float';
        if isempty(type_out), type_out = 'single'; end
    otherwise
        error([num2str(dtype) ' from get_head unrecognized data type'])
end

if framesize~=ni*nj*nbytes
    disp 'BAD HEADER!!! framesize does not match framewidth*frameheight*nbytes!'
    framesize = ni*nj*nbytes;
end

% DIFFERENTIAL CAMERA? (then each condition stores one reference image
% between the differential data)
includesrefframe = ( (filesize-lenh) > (framesize*nfr*nstim) );
if includesrefframe
    nfr2 = nfr+1;
else
    nfr2 = nfr;
end

% BLOCK BINNING PARAMETERS
if ischar(flag)
    switch flag
        case 'av'
            flag=nblocks;
        case ''
            flag=1;
        otherwise
            error('unknown flag ''%s''',flag)
    end
end
nbin=round(nblocks/flag);

% CONDITION FILTERING
if isempty(conds), conds = 1:nstim; end
ncond = length(conds);

% READING BLOCKS
if nblocks>1
    try fn_progress('reading block',nblocks), end %#ok<TRYNC>
end
% allocate
if doarray
    data = zeros(ni,nj,nfr,ncond,nbin,type_out);
else
    data = cell(nbin,ncond);
    for k=1:nbin
        for i=1:ncond
            data{k,i} = zeros(ni,nj,nfr,type_out);
        end
    end
end
% read
for k=1:nbin
    bin = floor((k-1)/nbin*nblocks)+1:floor(k/nbin*nblocks);
    sbin = length(bin);
    for j=1:sbin
        if nblocks>1
            try fn_progress(bin(j)), catch fprintf('reading block %i/%i\n',bin(j),nblocks), end
        end
        
        % open file
        [file, msg] = fopen(fname(bin(j),:), 'r', 'l'); %assume files are little-endian (PC or VAX)
        if (file == -1), error(msg); end;
        
        % read
        for i=1:ncond
            
            % % a = loadsum(deblank(),i*nfr2-nfr:i*nfr2-1,0,datatype,ni,nj,lenh);
            % do not call loadsum
            framestart = conds(i)*nfr2-nfr;
            offset = framestart * ni * nj * nbytes + lenh;
            fseek(file, offset, 'bof');
            a = fread(file,[ni*nj nfr],[datatype '=>' type_out]);
            
            if size(a,2)==nfr
                a = reshape(a,ni,nj,nfr);
            else
                disp 'BIG PROBLEM: number of frames does not fit headers'
                nfrok = size(a,2);
                a = reshape(a,ni,nj,nfrok);
                a = cat(3,a,repmat(mean(a,3),[1 1 nfr-nfrok]));
            end
            
            if includesrefframe
                framestart = (conds(i)-1)*nfr2;
                offset = framestart * ni * nj * nbytes + lenh;
                fseek(file, offset, 'bof');
                ref = fread(file,[ni nj],[datatype '=>' type_out]);
                a = repmat(ref*scalefactor,[1 1 nfr]) + a/gain - scalefactor*dc/gain;
            end
            if sbin==1
                if doarray
                    data(:,:,:,i,k) = a;
                else
                    data{k,i} = a;
                end
            else
                if doarray
                    data(:,:,:,i,k) = data(:,:,:,i,k) + a/sbin;
                else
                    data{k,i} = data{k,i} + a/sbin;
                end
            end
        end
        
        % close file
        fclose(file);
        
    end
end

