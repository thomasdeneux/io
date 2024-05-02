function [recordings vectors menupar xpar] = elphy_read(files)
% function [recordings vectors menupar xpar] = elphy_read(files)
%---
% warning: for now, reads only content from first channel

if nargin<1, files = brick.getfile('*.DAT'); end

files = cellstr(files);
nf = length(files);
if nf>1, error 'multiple files not implemented yet', end
[recordings vectors menupar xpar] = ReadOneFile(files{1});


function [recordings vectors menupar xpar] = ReadOneFile(fname)

if ~exist(fname,'file')
    error 'file does not exist'
end

% copy file to local directory for faster reading afterwards
tempDir = getenv('temp');
fbase = brick.fileparts(fname,'name'); % remember the original name of the file
if ~ismember(fname(1),'CD') && exist(tempDir,'dir')
    tempName = fullfile(tempDir,'Temp_LoadElphyBehavior.DAT');
    copyfile(fname,tempName)
    fname = tempName;
end

[recordings dates vectors menupar xpar epinfo] = Read(fname,fbase);
vectors.dates = dates;
if ~isempty(epinfo)
    F = fieldnames(epinfo);
    for k=1:length(F)
        f = F{k};
        vectors.(f) = [epinfo.(f)];
    end
end

% Fix some vectors: first get the appropriate length values
if isequal(xpar,struct)
    disp 'Typical parameters from Brice setup not found, finishing here.'
    return
end

nepmax = xpar.fix.MaxTrialNumber;
nstimperep = xpar.fix.NStimPerEp;
if iscell(recordings)
    nep=length(recordings);
else
    nep = size(recordings,2);
end

if isfield(vectors,'iCond') && length(vectors.iCond)~=nep
    error 'vector ''iCond'' is not of the expected length'
end

nstim = nep*nstimperep;
nstimmax = nepmax*nstimperep;

% when Elphy crashes in the middle of the experiment, vectors 'trecord',
% 'crecord' and 'lickrecord' are not saved; moreover, if the next
% experiment is stoped before the first episode completed, erroneous
% vectors will be saved in the file (instead of in the file of the next
% experiment; in this latter case the vectors are filled with zeros)
errmsg = '';
if ~isfield(vectors,'trecord')
    errmsg = 'no vector ''trecord'' was saved';
elseif any(vectors.trecord(1:size(recordings,2)*xpar.fix.NStimPerEp)==0)
    errmsg = 'vector ''trecord'' has zeros';
elseif length(vectors.trecord)==xpar.fix.MaxTrialNumber && xpar.fix.NStimPerEp>1
    warning 'Length of trecord is equal to MaxTrialNumber, whereas it should be multiplied by NStimPerEp. Not generating an error.'
elseif length(vectors.trecord)~=xpar.fix.MaxTrialNumber*xpar.fix.NStimPerEp
    errmsg = 'vector ''trecord'' is not of the expected length (STRANGE THAT IT DOES NOT HAVE ZEROS...)';
end
    
% Try to infer 'trecord' from 'iCond' or 'tprebuild' and 'crecord' from
% 'correct'
if ~isfield(vectors,'iCond') || (nstimperep>1 && ~isfield(vectors,'tprebuild'))
    fixmsg = '''trecord'', ''crecord'' and ''lickrecord'' will not be available';
    fixtrec = [];
    fixcrec = [];
elseif nstimperep>1
    % variables 'iCond' and 'correct' are saved only once per episode
    % (at the end of the episode), therefore they cannot be used to
    % recover vector 'trecord': we use 'tprebuild' instead
    if xpar.fix.DoRepeatCondition, error 'DoRepeatCondition=True and NStimPerEp>1 is not supposed to be possible', end
    if length(vectors.tprebuild)~=xpar.fix.MaxTrialNumber*xpar.fix.NStimPerEp, error 'vector ''tprebuild'' is not of the expected length', end
    fixmsg = '''trecord'' could be inferred from ''tprebuild''; ''crecord'' and ''lickrecord'' could not be recovered';
    fixtrec = vectors.tprebuild(1:nstim);
    fixcrec = [];
else
    % variables 'iCond' and 'correct' can be used, we check if
    % applicable that they match with tprebuild
    fixmsg = '''trecord'' and ''crecord'' could be inferred from ''iCond'' and ''correct''; ''lickrecord'' could not be recovered';
    fixtrec = vectors.iCond;
    fixcrec = vectors.correct;    
    if ~xpar.fix.DoRepeatCondition && isfield(vectors,'tprebuild') && ~isequal(fixtrec,vectors.tprebuild(1:nstim))
        error '''iCond'' and ''tprebuild'' do not match'
    end
end
if isempty(errmsg)
    vectors.trecord = vectors.trecord(1:nstim);
    vectors.crecord = vectors.crecord(1:nstim);
    vectors.lickrecord = vectors.lickrecord(1:nstim);
    % check that the saved trecord matches the inferred one if any
    if ~isempty(fixtrec) && ~isequal(fixtrec,vectors.trecord)
        if nstimperep>1
            error 'vectors ''trecord'' and ''tprebuild'' do not match'
        else
            error 'vectors ''trecord'' and ''iCond'' do not match'
        end
    end
else
    disp([fbase ': ' errmsg])
    disp([fbase ': ' fixmsg])
    vectors.trecord = fixtrec;
    vectors.crecord = fixcrec;
    vectors.lickrecord = zeros(1,nstim);
end



%------------- MAIN SUB-FUNCTIONS -----------------------------------------

function [recordings dates vectors menupar xpar epinfo] = Read(fname,fbase)

recordings = {};
vectors = struct;
menupar = [];
memonames = {};
xpar = struct;
epinfo = [];

% open file
fid=fopen(fname);
if (fid==-1)
    error 'unable to open file'
end

% skip headers
SkipHeader(fid)

% read: loop on blocks
krec = 0;
while true
    [ID SZ posend] = getBlockID(fid);
    %fprintf('block %s: size %i\n',ID,SZ)
    if SZ==-1, break, end
    
    switch ID
        case 'B_Ep'
            % preallocate
            krec = krec+1;
            if krec>length(recordings), recordings{2*krec} = []; end %#ok<AGROW>
            % read
            [recordings{krec} dates{krec}] = ReadEpisode(fid,SZ,posend); %#ok<AGROW>
        case 'Vector'
            [name value] = ReadVector(fid,SZ,posend);
            name = lower(name); name = strrep(name,'pg0.','');
            if isfield(vectors,name)
                disp([fbase ': several ''' name ''' vectors saved in the file! keeping only the first one'])
            else
                vectors.(name) = brick.row(value);
            end
            fseek(fid,posend,'bof');
        case 'DBrecord'
            if krec==0
                disp 'ignoring DBrecord before first episode'
            else
                p = ReadParameters(fid,SZ,posend);
                if isfield(p,'ProtocolType')
                    if isempty(menupar)
                        menupar = p;
                    else
                        disp 'problem! several set of parameters in file!'
                    end
                else
                    if isstruct(epinfo) && length(epinfo)==krec
                        fprintf('%s: variables saved several time for episode %i, keeping only the first time\n',fbase,krec)
                    else
                        if krec==1
                            epinfo = p;
                        else
                            epinfo(krec) = p; %#ok<AGROW>
                        end
                    end
                end
            end
        case 'Memo'
            [memo name] = ReadMemo(fid,SZ,posend);
            if ismember(name,memonames)
                disp([fbase ': several ''' name ''' memos saved in the file! keeping only the first one'])
            else
                memonames{end+1} = name; %#ok<AGROW>
                xpar = brick.structmerge(xpar,memo);
            end
        otherwise
            % this function is used only when trying to extract seldom
            % information from corrupted file
            if ~strcmp(ID,'DATA')
                xpar = ReadOtherField(fid,krec,ID,SZ,posend,xpar);
            end
    end
    
    fseek(fid,posend,'bof');
end

% close file
fclose(fid);

% rearrange recordings
if iscell(recordings{1})
    % several channels, with different number of samples
    disp 'channels do not all have same number of samples, returning a cell array rather than a matrix'
    recordings = recordings(1:krec);
elseif isvector(recordings{1})
    % only one channel
    try
        recordings = [recordings{:}];
    catch        
        disp 'episodes do not all have same number of samples, returning a cell array rather than a matrix'
        emptyrec = brick.isemptyc(recordings);
        recordings = recordings(1:find(~emptyrec,1,'last'));
    end
else
    % several channels
    recordings = permute(cat(3,recordings{:}),[1 3 2]);
end

function SkipHeader(fid)

% Elphy Header
len=fread(fid,1,'uchar');
ID=fread(fid,len,'uchar');
ID=char(ID');
dum=fread(fid,15-len,'uchar');
SZ=fread(fid,1,'int16');

if (len ~=12) || ( strcmp(ID ,'Dac2 objects'))
    error problem
end;


function [V date] = ReadEpisode(fid,SZ,PosMax) %#ok<INUSL>
% Read recordings

% read episode headers (sub-blocks of the block 'B_Ep')
while ftell(fid)<PosMax
    
    [ID SZ pos1] = subBlockID(fid);
    
    switch ID
        case 'Ep'
            nbvoie= fread(fid,1,'uchar');
            
            nbpt= fread(fid,1,'int32');
            
            tpData= fread(fid,1,'uchar');
            uX= readString(fid,10);
            
            Dxu = fread(fid,1,'double');
            x0u = fread(fid,1,'double');
            
            continu = fread(fid,1,'uchar');
            
            TagMode = fread(fid,1,'uchar');
            TagShift = fread(fid,1,'uchar');
            
            if SZ>36
                DxuSpk = fread(fid,1,'double');
                X0uSpk = fread(fid,1,'double');
                nbSpk = fread(fid,1,'int32');
                DyuSpk = fread(fid,1,'double');
                Y0uSpk = fread(fid,1,'double');
                unitXspk = readString(fid,10);
                unitYSpk = readString(fid,10);
                CyberTime = fread(fid,1,'double');
                PCtime = fread(fid,1,'int32');
            end;
            
            Ktype = ones(nbvoie,1) * 2 ;       %Default type = smallint
            
        case 'Adc'
            uY=[];
            Dyu=[];
            Y0u=[];
            for ii = 1 : nbvoie
                uY = strvcat(uY,readString(fid,10));
                Dyu = cat(1,Dyu,fread(fid,1,'double'));
                Y0u = cat(1,Y0u,fread(fid,1,'double'));
            end;
            
        case 'Ksamp'
            Ksamp = [];
            for ii = 1 : nbvoie
                Ksamp = cat(1,Ksamp,fread(fid,1,'uint16'));
            end;
            
        case 'Ktype'
            Ktype=[];
            for ii = 1 : nbvoie
                Ktype = cat(1,Ktype,fread(fid,1,'uchar'));
                
            end
    end
    
    fseek(fid,pos1,'bof');
end
fseek(fid,PosMax,'bof');

% read data ('RDATA' block)
SZ = findBlock(fid,'RDATA');
if SZ==-1
    error 'no episode data'
end
RdataHsize = fread(fid,1,'int16');    %Rdata Header Size
% fseek(fid,RdataHsize-2,'cof');
curpos = ftell(fid);
fread(fid,1); 
x=fread(fid,1,'uint64');
date = datestr(x*2^-37-33286456);
fseek(fid,curpos+RdataHsize-2,'bof'); % let's be on the safe size for repositionning
SZ = SZ-RdataHsize;

if SZ==0
    disp 'ecountered empty recording'
    V = [];
    return
end
nbSamp = floor(SZ/2);

[AgSampleCount,ppcm0] = GetAgSampleCount(Ksamp);
ChanMask = GetMask(Ksamp);
SamplePerChan = GetSamplePerChan(nbSamp,AgSampleCount,ppcm0,Ksamp,Ktype,ChanMask);

N = SamplePerChan;
Dy0 = Dyu;
Y00 = Y0u;

if (TagMode == 1)
    Dy0 = Dy0/power(2,TagShift);
end;
Dy0(Ktype==5) = 1;

if any(diff(Ktype==5))
    if nbvoie>1, error 'not implemented yet', end
    % different types of numbers, go the slow way
    % BEGIN OLD CODE
    V = zeros(N,1);
    k =0;
    for i=1 : nbSamp
        im = mod((i-1),AgSampleCount)+1;
        if Ktype(ChanMask(im))==5
            w = fread(fid,1,'float32');
        else
            w = fread(fid,1,'int16');
        end;
        if ChanMask(im) == NumChan
            k = k+1;
            V(k) = w*Dy0 + Y00;
        end;
    end;
    % END OLD CODE
else
    % go the fast way: load everything at once!
    if Ktype(1)==5
        w = fread(fid,[AgSampleCount nbSamp/AgSampleCount],'float32');
    else
        w = fread(fid,[AgSampleCount nbSamp/AgSampleCount],'int16');
    end;
    V = cell(1,nbvoie);
    for i=1:nbvoie
        V{i} = brick.column(w(ChanMask==i,:))*Dy0(i) + Y00(i);
    end
    if ~any(diff(N)) % all channels have same number of records
        V = [V{:}];
    end
end


function [name V] = ReadVector(fid,SZ,PosMax) %#ok<INUSL>

name = [];

% Read header (from sub-blocks of the 'Vector' block)
while ftell(fid)<PosMax
    
    [ID SZ pos1] = subBlockID(fid); 
    switch ID
        case 'IDENT1'
            name = fread(fid,[1 SZ],'*char');
        case 'OBJINF'
            tpNum = fread(fid,1,'uchar');

            imin = fread(fid,1,'int32');
            imax = fread(fid,1,'int32');
            jmin = fread(fid,1,'int32'); % not used for Tvector
            jmax = fread(fid,1,'int32'); % not used for Tvector

            x0u = fread(fid,1,'float64');
            dxu = fread(fid,1,'float64');
            y0u = fread(fid,1,'float64');
            dyu = fread(fid,1,'float64');
    end
    fseek(fid,pos1,'bof');
    
end
fseek(fid,PosMax,'bof');

% Read data ('DATA' block)
SZ = findBlock(fid,'DATA');
if SZ==-1
    error 'no data for vector'
end

% (we must skip one byte)
fseek(fid,1,'cof');

% (then read data)
format = Nformat(tpNum);
V = fread(fid, imax-imin+1, format);

% Apply scaling params
V = V*dyu+y0u;


function parameters = ReadParameters(fid,SZ,PosMax) %#ok<INUSL>
% Read header (from sub-blocks of the 'Vector' block)

while ftell(fid)<PosMax
    [ID SZ pos1] = subBlockID(fid); 
    switch ID
        case 'ST'
            str = fread(fid,[1 SZ],'*char');
            sep=[-1 strfind(str,char([13 10]))]; 
            npar = length(sep)-1; 
            names = cell(1,npar);
            for i=1:npar, names{i}=str(sep(i)+2:sep(i+1)-1); end
        case 'BUF'
            values = cell(1,npar);
            for i=1:npar
                type = fread(fid,1,'uchar');
                switch type
                    case 1
                        values{i} = fread(fid,1,'uchar=>logical');
                    case 2
                        values{i} = fread(fid,1,'int64');
                    case 3
                        %values{i} = extended2double(fread(fid,10,'uchar'));
                        x = fread(fid,1,'uint64')*2^-63;
                        s_e = fread(fid,1,'int16'); % first bit is the sign, the rest is the exponent coded as an int15
                        s = sign(s_e);
                        e = abs(s_e)-(2^14-1);
                        values{i} = sign(s_e)*x*2^e;
                    case 4
                        len = fread(fid,1,'uint32');
                        if len, values{i} = fread(fid,len,'*char')'; else values{i}=''; end
                    case 7
                        len = fread(fid,1,'uchar');
                        values{i} = fread(fid,[1 len],'*char')';
                    otherwise
                        disp(sprintf('unknown type number %i, skipping parameter %s',type,names{i})) %#ok<DSPS>
                end
            end
    end
    fseek(fid,pos1,'bof');
    
end
fseek(fid,PosMax,'bof');

C = [names; values];
parameters = struct(C{:});


function [xpar memoname] = ReadMemo(fid,SZ,PosMax) %#ok<INUSL>

while ftell(fid)<PosMax
    [ID SZ pos1] = subBlockID(fid); 
    switch ID
        case 'IDENT1'
            memoname = fread(fid,[1 SZ],'*char');
        case 'ST'
            memo = fread(fid,SZ,'*char');
        otherwise
            %             disp([ID ' - ' num2str(SZ)])
    end
    fseek(fid,pos1,'bof');
    
end
fseek(fid,PosMax,'bof');

memo = memo';

% Convert Memo into structured information
switch memoname
    case 'PG0.PPAR2'
        s = struct;
        isunique = struct;
        lines = brick.strcut(memo,[10 13]);
        for i=1:length(lines)
            items = brick.strcut(lines{i},';');
            name = items{1};
            isunique.(name) = (length(items)<3 || isempty(items{3}));
            if isunique.(name)
                if length(items)<2, values = {''}; else values = items(2); end
            else
                values = items(2:end);
            end
            for j=1:length(values)
                if strcmp(values{j},'FALSE')
                    values{j} = false;
                elseif strcmp(values{j},'TRUE')
                    values{j} = true;
                else
                    val = str2double(values{j});
                    if ~isnan(val), values{j} = val; end
                end
            end
            s.(name) = values;
        end
        isunique = cell2mat(struct2cell(isunique));
        names = fieldnames(s);
        values = struct2cell(s);
        if all(isunique)
            xpar = struct( ...
                'fix',      cell2struct([values{isunique}]',names(isunique)), ...
                'table',    struct);
        else
            xpar = struct( ...
                'fix',      cell2struct([values{isunique}]',names(isunique)), ...
                'table',    cell2struct(cat(1,values{~isunique}),names(~isunique))');
        end
    case 'PG0.IMAGELIST'
        xpar = struct('imagelist',{brick.strcut(memo,[10 13])});
    case 'PG0.SOUNDLIST'
        xpar = struct('soundlist',{brick.strcut(memo,[10 13])});
    otherwise
        fprintf('don''t know how to read Memo %s',memoname)
        xpar = struct(memoname,memo);
end


%------------- ROUTINE SUB-FUNCTIONS --------------------------------------


function [SZ posend] = findBlock(fid,Name)
ID=0;
SZ=0;
N=0;

while N==0 && ~feof(fid)
    SZ=fread(fid,1,'int32');
    if (length(SZ)==0 )
        break;
    end;
    
    len=fread(fid,1,'uchar');
    if (length(len)==0 )
        break;
    end;
    
    ID=fread(fid,len,'uchar');
    ID=char(ID');
    
    if ( strcmp(ID , Name))
        N=N+1;
    end;
    SZ=SZ-len-5;
    if N==0
        fseek(fid,SZ,'cof');
    end;
end;
if N==0
    SZ=-1;
    %disp 'Block not found'
end

posend = ftell(fid)+SZ;


function [ID SZ posend] = getBlockID(fid)

SZ=fread(fid,1,'int32');
if isempty(SZ), SZ=-1; ID=[]; posend=[]; return, end
    
len=fread(fid,1,'uchar');
if isempty(len), error problem, end
    
ID=fread(fid,[1 len],'*char');
    
SZ=SZ-len-5;
posend = ftell(fid)+SZ;


function [ID, SZ, posend] = subBlockID(fid)
% read sub-block ID and size

len = fread(fid,1,'uchar');
if isempty(len)
    error problem
end;

ID=fread(fid,[1 len],'*char');
SZ=fread(fid,1,'uint16');
if (SZ==65535)
    SZ=fread(fid,1,'int32');
end;
posend = ftell(fid)+SZ;


function x = readString(fid,N)
% if N is specified, the cursor must be located at currentpos+1+N after the
% function returns
if nargin<2, N=0; end

len = fread(fid,1,'uchar');
if N, len = min(len,N); end

x = fread(fid,[1 len],'*char');

if N, fseek(fid,N-len,'cof'); end


function [x,ppcm0] = GetAgSampleCount(KS)
%Number of samples in aggregate

N = length(KS);
ppcm0 = 1;
for i =1 : N
    if KS(i)>0
        ppcm0 = lcm(ppcm0,KS(i));
    end;
end;

x = 0;
for i=1 : N
    if KS(i)>0
        x = x+ppcm0 / KS(i);
    end;
end;


function x = GetMask(KS)
% build Channel Mask = array of channel numbers for each sample of aggregate
% KS = array of DownSampling factors
AgC = GetAgSampleCount(KS);
Nvoie=length(KS);

x = zeros(AgC,1);
i=0;
k=1;
while k<=AgC
    for j=1 : Nvoie
        if (KS(j)>0) && (mod(i,KS(j))==0)
            x(k)=j;
            k=k+1;
            if k>AgC
                break;
            end;
        end;
    end;
    i=i+1;
end;


function x = GetAgSize(KS,KT)
tpSize =[1 1 2 2 4 4 6 8 10 8 16 20 0];
x=0;
for i= 1:length(KS)
    x = x+TpSize(KT(KS(i)));
end;


function x = GetSamplePerChan(nb,AgSz,ppcm0,KS,Ktype,chanMask)
% nb est la taille du bloc de donn�es
% AgSz est la taille de l'agr�gat
% KS contient les DownSampling factors

tpSize =[1 1 2 2 4 4 6 8 10 8 16 20 0];

x = zeros(length(KS),1);
nvoie = length(KS);

nbAg = floor(nb/AgSz); %        {nombre d'agr�gats complets }

for i=1 : nvoie
    if KS(i)>0
        x(i) = nbAg*ppcm0/KS(i);
    else
        x(i) = 0;
    end;
end;

rest= mod(nb,AgSz);        %   {agr�gat incomplet }
it=0; %{ taille }
j=1;  %{ indice dans l'Ag }
while it<rest
    vv = chanMask(j);
    if (vv>0) && (vv <=nvoie)
        x(vv) = x(vv)+1;
        it = it + tpSize(Ktype(vv));
        j = j+1;
    end;
end;


function Nf = Nformat(N)
switch N
    case 0
        Nf ='uint8';
    case 1
        Nf ='int8';
    case 2
        Nf ='int16';
    case 3
        Nf ='uint16';
    case 4
        Nf ='int32';
    case 5
        Nf ='float32';
    case 7
        Nf ='float64';
    case 9
        Nf ='float32';
    case 10
        Nf ='float64';
    otherwise
        Nf = '';
end;


function x=extended2double(e)

bits = dec2bin(e,8);
bits = reshape(flipud(bits)',1,80);

sign = (bits(1)=='1');
exponent = bin2dec((bits(2:16)));

% normalizeCorrection = (bits(17)=='1');
% m = (bits(18:80));
% mantissa = bin2dec(m(1:52)) * 2^11 + bin2dec(m(53:63)); % note that precision in m(53:63) is actually lost
% x = (-1)^sign * (normalizeCorrection + mantissa*2^-63) * 2^(exponent-16383);

u = bin2dec(bits(17:68)) * 2^(12-63);
x = (-1)^sign * u * 2^(exponent-16383);

%------------- MANUAL EXTRACTION OF SPECIFIC DATA -------------------------


function xpar = ReadOtherField(fid,krec,ID0,SZ,PosMax,xpar) %#ok<INUSL>

% % which image file?
% if ~strcmpi(ID0,'VSbitmap'), return, end
% 
% while ftell(fid)<PosMax
%     [ID SZ pos1] = subBlockID(fid);   %#ok<ASGLU>
%     if strcmp(ID,'StFile')
%         str = fread(fid,[1 SZ],'*char');
%         xpar.stfile{krec} = brick.fileparts(str,'name');
%     end
%     fseek(fid,pos1,'bof');
%     
% end
% fseek(fid,PosMax,'bof');

return
% visual object
if ~brick.ismemberstr(ID0,{'Bar' 'LGrating' 'disk'}), return, end

while ftell(fid)<PosMax
    [ID SZ pos1] = subBlockID(fid);   %#ok<ASGLU>
    p = ftell(fid);
    switch ID0
        case 'Bar'
            if strcmp(ID,'IDENT1')
                name = fread(fid,[1 SZ],'*char');
                if ~strcmp(name,'PG0.BAR'), return, end
            elseif strcmp(ID,'DEG')
                fseek(fid,22,'cof');
                lum = fread(fid,1,'single');
                xpar.barlum(krec) = lum;
            end
        case 'LGrating'
            if strcmp(ID,'ORIENT')
                oris = [0 45 90 135 180 225 270 315];
                ori = fread(fid,1,'single');
                xpar.gtype(krec) = find(ori==oris,1);
            end
        case 'disk'
            if strcmp(ID,'DEG')
                degs = [2 16];
                deg = fread(fid,4,'single');
                xpar.dtype(krec) = find(deg(3)==degs,1);
            end
    end
    fseek(fid,pos1,'bof');
    
end
fseek(fid,PosMax,'bof');


