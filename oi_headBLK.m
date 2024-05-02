function p = oi_headBLK(filename)
% function p = oi_headBLK(filename)
%---
% returns structure with all parameter values of a block file

if nargin<1 
    try
        filename = fn_getfile('*.BLK','Select block file');
    catch
        [filename pathname] = uigetfile('*.BLK','Select block file','MultiSelect','on');
        filename = fullfile(pathname,filename);
        filename = char(filename);
    end
end
    
fid=fopen(filename,'r','l');
if fid == -1,
	error(['Invalid filename: ',filename])
end
%%%%%%%%%%%%%%%%%% BEGIN DEFINITIONS OF VARIABLES %%%%%%%%%%%%%%%%%%

% note on data type:
% long = int32
% ushort = uint16
% float = float32
% char = (uint8=>char), i.e. read a uint8 and convert to char; reading
% directly a char makes a problem because
% for example, if byte is read as a uint8 as 220, (uint8=>char) reading
% format results in char(220), but char format results in 65533, and
% (char=>char) format results in char(65533)

% DATA INTEGRITY 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.filesize              =fread(fid,1,'int32');
p.checksum_header			=fread(fid,1,'int32');
	% beginning with the lLen Header field

p.checksum_data			=fread(fid,1,'int32');

% COMMON TO ALL DATA FILES 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.lenheader			=fread(fid,1,'int32');
p.versionid			=fread(fid,1,'float32');
p.filetype			=fread(fid,1,'int32');
	% RAWBLOCK_FILE          (11)
	% DCBLOCK_FILE           (12)
	% SUM_FILE               (13)

p.filesubtype			=fread(fid,1,'int32');
	% FROM_VDAQ              (11)
	% FROM_ORA               (12)
	% FROM_DYEDAQ            (13)

p.datatype			=fread(fid,1,'int32');
	% DAT_UCHAR     (11)
	% DAT_USHORT    (12)
	% DAT_LONG      (13)
	% DAT_FLOAT     (14)

p.sizeof				=fread(fid,1,'int32');
	% e.g. sizeof(long), sizeof(float)

p.framewidth			=fread(fid,1,'int32');
p.frameheight			=fread(fid,1,'int32');
p.nframesperstim			=fread(fid,1,'int32');
p.nstimuli			=fread(fid,1,'int32');
p.initialxbinfactor		=fread(fid,1,'int32');
	% from data acquisition
p.initialybinfactor		=fread(fid,1,'int32');
	% from data acquisition

p.xbinfactor			=fread(fid,1,'int32');
	% this file
p.ybinfactor			=fread(fid,1,'int32');
	% this file

p.username			=fread(fid,32,'uint8=>char')';
p.recordingdate			=fread(fid,16,'uint8=>char')';
p.x1roi				=fread(fid,1,'int32');
p.y1roi				=fread(fid,1,'int32');
p.x2roi				=fread(fid,1,'int32');
p.y2roi				=fread(fid,1,'int32');

% LOCATE DATA AND REF FRAMES 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.stimoffs			=fread(fid,1,'int32');
p.stimsize			=fread(fid,1,'int32');
p.frameoffs			=fread(fid,1,'int32');
p.framesize			=fread(fid,1,'int32');
p.refoffs				=fread(fid,1,'int32');
p.refsize				=fread(fid,1,'int32');
p.refwidth			=fread(fid,1,'int32');
p.refheight			=fread(fid,1,'int32');

% Common to data files that have undergone some form of
% "compression" or "summing"
% i.e. The data in the current file may be the result of
%      having summed blocks 'a'-'f', frames 1-7
p.whichblocks			=fread(fid,16,'uint16');
p.whichframes			=fread(fid,16,'uint16');

% DATA ANALYSIS			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.loclip				=fread(fid,1,'int32');
p.hiclip				=fread(fid,1,'int32');
p.lopass				=fread(fid,1,'int32');
p.hipass				=fread(fid,1,'int32');
p.operationsperformed		=setstr(fread(fid,64,'uint8=>char'))';

% ORA-SPECIFIC			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.magnification			=fread(fid,1,'float32');
p.gain				=fread(fid,1,'uint16');
p.wavelength			=fread(fid,1,'uint16');
p.exposuretime			=fread(fid,1,'int32');
p.nrepetitions			=fread(fid,1,'int32');
p.acquisitiondelay		=fread(fid,1,'int32');
p.interstiminterval		=fread(fid,1,'int32');
p.creationdate			=setstr(fread(fid,16,'uint8=>char'))';
p.datafilename			=setstr(fread(fid,64,'uint8=>char'))';
p.orareserved			=setstr(fread(fid,256,'uint8=>char'))';


if p.filesubtype == 13,   %it's dyedaq file
  
  %  OIHEADER.H
  %  last revised 4.5.97 by Chaipi Wijnbergen for DyeDaq
  %  
  %  DyeDaq-specific
  p.includesrefframe =fread(fid,1, 'int32');     % 0 or 1
  p.temp =fread(fid,128, 'uint8=>char');
  p.listofstimuli=temp(1:max(find(temp~=0)))';  % up to first non-zero stimulus
  p.ntrials =fread(fid,1, 'int32');
  p.scalefactor =fread(fid,1, 'int32');          % bin * trials
  p.cameragain =fread(fid,1, 'short');         % shcameragain        1,   2,   5,  10
  p.ampgain =fread(fid,1, 'short');            % amp gain            1,   4,  10,  16,
                                             %                    40,  64, 100, 160,
                                             %                    400,1000
  p.samplingrate =fread(fid,1, 'short');       % sampling rate (1/x)
                                             %                     1,   2,   4,   8,
                                             %                     16,  32,  64, 128,
                                             %                     256, 512,1024,2048
  p.average =fread(fid,1, 'short');            % average             1,   2,   4,   8,
                                             %                    16,  32,  64, 128
  p.exposuretime =fread(fid,1, 'short');       % exposure time       1,   2,   4,   8,
                                             %                    16,  32,  64, 128,
                                             %                    256, 512,1024,2048
  p.samplingaverage =fread(fid,1, 'short');    % sampling average    1,   2,   4,   8,
                                             %                    16,  32,  64, 128
  p.presentaverage =fread(fid,1, 'short');
  p.framesperstim =fread(fid,1, 'short');
  p.trialsperblock =fread(fid,1, 'short');
  p.sizeofanalogbufferinframes =fread(fid,1, 'short');
  p.cameratrials=fread(fid,1, 'short');
  p.filler =setstr(fread(fid,106, 'uint8=>char'))';
  
  p.dyedaqreserved =setstr(fread(fid,256, 'uint8=>char'))';
else   % it's not dyedaq specific

  % VDAQ-SPECIFIC			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p.includesrefframe		=fread(fid,1,'int32');
  p.listofstimuli			=setstr(fread(fid,256,'uint8=>char'))';
  p.nvideoframesperdataframe	=fread(fid,1,'int32');
  p.ntrials				=fread(fid,1,'int32');
  p.scalefactor			=fread(fid,1,'int32');
  p.meanampgain			=fread(fid,1,'float32');
  p.meanampdc			=fread(fid,1,'float32');
  p.vdaqreserved			=setstr(fread(fid,256,'uint8=>char'))';
end    % end of VDAQ specific

% USER-DEFINED			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p.user				=setstr(fread(fid,256,'uint8=>char'))';

% COMMENT			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.comment				=setstr(fread(fid,256,'uint8=>char'))';
p.refscalefactor =fread(fid,1, 'int32');          % bin * trials for reference

%%%%%%%%%%%%%%%%%%  END DEFINITIONS OF VARIABLES  %%%%%%%%%%%%%%%%%%

% Note that doing ftell(fid) here results in p.lenheader + 4!!, which means
% that we read one extra int32!!! But which one...?

fseek(fid,0,1);            % go to EOF
p.actuallength=ftell(fid);   % report where EOF is in bytes

fclose(fid);

