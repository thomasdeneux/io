function [desc data sfr times stims stimunext] = tpm_2plsmread(fname,verbose,frnum)
% function [desc data times stims stimunext] = tpm_2plsmread(fname,verbose,frnum)


if ~exist('verbose','var')
    verbose = 0;
end

fid = fopen(fname,'r','b');
% 'b' stands for big-endian byte ordering
% fid is the .2plsm file containing all the data

fseek(fid,0,1);
endoffile = ftell(fid);
fseek(fid,0,-1);

cnt = 1; acnt = 1;
done = 0;

while ~done %&& acnt<2
    
    bnel = fread(fid,1,'uint32');
    % bnel is the number of channels saved in SaveFrame.vi (14)
    dat = unflatten2string(fid,'char',verbose);
    tim = unflatten2string(fid,'char',verbose);
    nel = fread(fid,1,'uint32');
    gg = fread(fid,16,'uchar');
    nbytes = fread(fid,1,'uint32');
    nel = fread(fid,1,'uint32');
    
    %
    % parnames/-vals: names and values of frames parameters
    parnames = [];
    for j=1:nel
        parnames{j} = unflatten2string(fid,'char',verbose);
    end
    nel = fread(fid,1,'uint32');
    parvals = unflatten2string(fid,'double',verbose);
    % so de = desc contains the technical info about the frames
    de.AOrate = 0;
    for j=1:length(parvals)
        if cnt==1
            d =setfield(de,rmwhitespace(parnames{j}),parvals(j));
            de = d;
        else
            de=setfield(de,rmwhitespace(parnames{j}),parvals(j));
        end
    end
    de.frnum = acnt;
    %
    
    
    %
    % mdat = the 2ph imaging data
    mdat = unflatten2string(fid,'2darray=>single',verbose);
    %
    
    
    %
    % Here is the 3rd data channel(?)
    if bnel>6
        vdata = unflatten2string(fid,'1darray',verbose);
    end
    
    if bnel>7
        fst = unflatten2string(fid,'char');
    end
    
    if bnel>8
        nbytes = fread(fid,1,'uint32');
        rt = unflatten2string(fid,'int32');
    end
    
    if bnel>9
        % this channel, 10/14 contains the electrophy data
        vdat{cnt} = unflatten2string(fid,'1darray');
    end
    
    % channels 11 to 14 contain: comments (11), stimuli (12)
    % times of frames and stimuli (?, 13 & 14)
    
    % comments
    if bnel>10
        vvar1 = unflatten2string(fid,'char');
    end
    
    if bnel>11
        %      vvar2 = unflatten2string(fid,'char');
        strlen = fread(fid,1,'uint32');
        nel = fread(fid,1,'uint32');
        if nel > 0
            if nel == 1
                vvar2 = unflatten2string(fid,'char',verbose);
            else
                vvar2=cell(1,nel);
                for j=1:nel
                    vvar2{j} = unflatten2string(fid,'char',verbose);
                end
            end
        else
            vvar2 = 0;
        end
    end
    
    if bnel>12
        strlen = fread(fid,1,'uint32');
        nel = fread(fid,1,'uint32');
        if nel > 0
            vvar3 = zeros(1,nel);
            for j=1:nel
                vvar3(j) = fread(fid, 1, 'uint32');
            end
        else
            vvar3 = 0;
        end
    end
    
    if bnel>13
        sizestr = fread(fid, 1, 'uint32');
        vvar4 = fread(fid, 1, 'uint32');
    end
    
    %   for cel = 11:bnel
    %     unflatten2string(fid,'unknown');
    %   end
    
    
    %
    
    
    %
    % filling output data variables
    tdata =0;
    if exist('frnum','var')
        if sum(acnt==frnum)>0
            tdata = 1;
            if acnt == frnum(end)
                done = 1;
            end
        end
    end
    if tdata || ~exist('frnum','var')
        % First allocate some data
        if cnt == 1
            if ~exist('frnum','var')
                nrecs = floor(endoffile/ftell(fid));
            else
                nrecs = length(frnum);
            end
            %disp([nrecs parvals(2) parvals(3)+1])
            fn_progress('loaded frame',nrecs);
            % nrecs is the number of frames, the 2 others, images dimensions
            nchan = size(mdat,1); % practically nchan=2: green and red
            % prepare for following preprocessings
            % - remove few columns on the left (image is often strongly deformed on the left)
            % - remove columns on the left (flyback)
            % - remove last line
            % - invert sign (strange that we need to do that)
            nx0 = parvals(2);
            ny0 = parvals(3)+1;
            xstart = round(nx0/20)+2;
            xend   = round(nx0*3/4)+1;
            nx = xend-xstart+1;
            ny = ny0-1;
            alldata = zeros(2,nx,ny,nrecs);
        end
        
        desc(cnt) = de; %#ok<*AGROW>
        fstr{cnt} = fst; %#ok<*NASGU>
        vdat{cnt} = vdata;
        rtimes{cnt} = rt;
        
        var{1,cnt} = vvar1;
        var{2,cnt} = vvar2;
        var{3,cnt} = vvar3;
        var{4,cnt} = vvar4;
        
        % preprocessing - TODO: read directly the data cut!!
        frame = reshape(mdat,2,nx0,ny0);
        frame = -frame(:,xstart:xend,1:ny);
        % store data
        alldata(:,:,:,cnt) = frame;
        
        fn_progress(cnt)
        
        cnt = cnt + 1;
        
    end
    if ftell(fid)==endoffile
        done = 1;
    else
        acnt = acnt + 1; % acnt goes up to reach the number of frames
    end
    %
    
end % while ~done
fclose(fid);


data = reshape(alldata(1,:),[nx ny nrecs]);
sfr  = reshape(alldata(2,:),[nx ny nrecs]);

% stimulation parameters
[times,stims,stimunext] = timeextract(var);

%---
function str = rmwhitespace(str)

whitespace = [9 10 11 12 13 32];
ind = ismember(str,whitespace);
str(ind) = [];

%---
function [times,stims,stimunext] = timeextract(var)
ss = size(var,2);
stims1 = [];
stims2 = [];
times = zeros(1,ss);
stimunext = zeros(1,ss);
% stimunext is numcond if there is a stimulus in next frame and 0 else
% used in movie display
time0 = var{4,1};
for i=1:ss
    times(i) = var{4,i};
    if var{3,i} == 0
    else
        stims1 = cat(2,stims1,var{3,i});
        if iscell(var{2,i})
            len = length(var{2,i});
            L = zeros(1,len);
            for k=1:len
                L(k) = findnb(var{2,i}{1,k});
            end
        else
            L = findnb(var{2,i});
        end
        stims2 = cat(2,stims2,L);
        if i>1
            stimunext(i-1) = L(1);
        end
    end
end
times = (times-time0)./1000;
stims1 = (stims1-time0)./1000;
stims = cat(1,stims1,stims2);

%---
function nb = findnb(st)
j = 0;
while st(j+1) ~= ' '
    j = j + 1;
end
nb = str2double(st(1:j));


