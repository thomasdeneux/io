function info = oi_infoBLKname(files)
% function info = oi_infoBLKname(files)
%---
% get information on time, condition, experiment and block, validity from
% BLK file name

if nargin==0, help oi_infoBLKname, return, end

files = cellstr(files);

% Analyze first file name to auto-detec fields
file = fn_fileparts(files{1},'base');

[idxs tok]=regexp(fliplr(file),'(\d+B)*(\d+E)*_*(\d{6})*(\d*C)*(T|F)*','tokenExtents','tokens');
if isempty(idxs)
    disp 'could not interpret file name'
    info = struct;
end
idxs = length(file)-(fliplr(idxs{1})-1);
tok = tok{1};

ok.exp = ~isempty(tok{2});
if ok.exp
    idx.exp = idxs(2,1)+1:idxs(2,2);
end
ok.block = ~isempty(tok{1});
if ok.block
    idx.block = idxs(1,1)+1:idxs(1,2);
end
ok.valid = ~isempty(tok{5});
if ok.valid
    idx.valid = idxs(5,1);
end
ok.cond = ~isempty(tok{4});
if ok.cond
    idx.cond = idxs(4,1)+1:idxs(4,2);
end
switch length(tok{3})
    case 0
        ok.date = false;
        ok.time = false;
    case 12
        ok.date = true;
        ok.time = true;
        idx.date = idxs(3,1)+(0:5);
        idx.time = idxs(3,1)+(6:11);
    case 6
        ex = str2num(reshape(fliplr(tok{3}),2,3)')'; %#ok<ST2NM>
        if any(ex>[24 60 60])
            ok.date = true;
            ok.time = false;
            idx.date = idxs(3,1)+(0:5);
        elseif any(ex(1:2)>[31 12])
            ok.date = false;
            ok.time = true;
            idx.time = idxs(3,1)+(0:5);
        else
            % sorry but we can't know if what we see is a date or a time!
            ok.date = false;
            ok.time = false;
        end
end

% Read every file
info = struct;
info(length(files)) = struct;
for i=1:length(files)
    file = fn_fileparts(files{i},'base');
    
    if ok.exp
        info(i).exp = str2double(file(idx.exp));
    end
    if ok.block
        info(i).block = str2double(file(idx.block));
    end
    if ok.valid
        info(i).valid = (file(idx.valid)=='T');
    end
    if ok.cond
        info(i).cond = str2double(file(idx.cond));
    end
    if ok.date
        date = file(idx.date);
        date = [date(1:2) '/' date(3:4) '/' date(5:6)];
        info(i).date = date;
    else
        date = '';
    end
    if ok.time && ok.date
        time = file(idx.time);
        time = [time(1:2) ':' time(3:4) ':' time(5:6)]; 
        info(i).time = time;
        %info(i).timenum = datenum([date ' ' time],'dd/mm/yy HH:MM:SS');
    end   
end