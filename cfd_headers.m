function headers = cfd_headers(filename)
% function header = cfd_headers(filename)

% Input
if nargin<1
    filename = fn_getfile('*.CFD');
end

% Open file
fid = fopen(deblank(filename),'r');
if fid==-1, error('could not open file ''%s''',filename), end

headers = struct('CFH',[],'ACV',[],'ASV',[],'WFM',[]);

% First part - CFH - CfNT Header Structure - 256 bytes (was 96)
% first number in comments is the starting byte within the structure
% fprintf('CFH: %i\n',ftell(fid))

headers.CFH.hdrsize     = []; 
headers.CFH.version     = fread(fid,1,'int16');         % 0 - program version (1.582 = 4582) (0)
headers.CFH.name        = fread(fid,14,'uchar=>char')';  % 2 - data file name (2)
headers.CFH.user        = fread(fid,16,'uchar=>char')';  % 16 - user name (16)
headers.CFH.time        = fread(fid,16,'uchar=>char')';  % 32 - time of acquisition (32)
headers.CFH.date        = fread(fid,16,'uchar=>char')';  % 48 - date of acquisition (48)
headers.CFH.start       = fread(fid,1,'int32');         % 64 - time at start of grab (64)
headers.CFH.stop        = fread(fid,1,'int32');         % 68 - time at end of grab (68)
headers.CFH.xpos        = [];
headers.CFH.ypos        = [];
headers.CFH.zpos        = fread(fid,1,'int32');         % 72 - z starting position in 1000's of microns (72)
headers.CFH.hdrsize     = fread(fid,1,'int32');         % 76 - length of header static + dynamic (76)
headers.CFH.xpos        = fread(fid,1,'int32');         % 80 - x starting position in nm (80)
headers.CFH.ypos        = fread(fid,1,'int32');         % 84 - y starting position in nm 
% record last xyz position too - MUM 020704
headers.CFH.endPosX     = fread(fid,1,'int32');         % 88 - last x position
headers.CFH.endPosY     = fread(fid,1,'int32');         % 92 - last y position
headers.CFH.endPosZ     = fread(fid,1,'int32');         % 96 - last z position
% Intensity parameters shifted 12 byte behind - MUM 020704
headers.CFH.lAttenuationDepth = fread(fid,1,'single');  % 100 - For recalculation of used intensity - new 011107 MUM
headers.CFH.bIntensityDevUsed = fread(fid,1,'int32');   % 104 - Has intensity controlling devive been used (102)
headers.CFH.startIntensityInPercent = fread(fid,1,'single');    % 108 - necessary to recalculate intensity used for single images - MUM 030108
fread(fid,36,'int32');        % headers.CFH.rsvd        = 112 - reserved bytes: 256 total - 112 used = 144 bytes = 36 longs each with 4 byte - MUM 030108
    
% Second part - acv
% fprintf('ACV: %i\n',ftell(fid))

headers.ACV.buflen      = fread(fid,1,'int32');         %  256 - 1 actual D/A buffer length, as calc'd    - type changed long -> DWORD - MUM 040512
headers.ACV.pixbuflen   = fread(fid,1,'int32');         %  260 - 2 actual pixel (dout) buffer length        - type changed long -> DWORD - MUM 040512
headers.ACV.flags.pol         = fread(fid,1,'ubit1');   % +bit 0        Bipolar/Unipolar, is now always bipolar - MUM 030916
headers.ACV.flags.clksrc      = fread(fid,1,'ubit1');   % +bit 1        Clock source for both D/A and Dout
headers.ACV.flags.waveform    = fread(fid,1,'ubit2');   % +bit 23        waveform type with 1 bit to grow
headers.ACV.flags.trigsrc     = fread(fid,1,'ubit1');   % +bit 4        Trigger source for both D/A and Dout
headers.ACV.flags.trigpol     = fread(fid,1,'ubit1');   % +bit 5        Trigger polarity when external
headers.ACV.flags.trigmod     = fread(fid,1,'ubit2');   % +bit 67        Trigger mode
headers.ACV.flags.lineScan    = fread(fid,1,'ubit1');   % +bit 8        Line Scan Mode
headers.ACV.flags.zt          = fread(fid,1,'ubit1');   % +bit 9        z or time series, stage or a dummy is active, replaces fTacqTimeBase - MUM 030708
headers.ACV.flags.bReturnToStart = fread(fid,1,'ubit1');% +bit a
fread(fid,1,'ubit1');   % headers.ACV.flags.fWrtEachImg = +bit b        write each image, now unused - MUM 040304
headers.ACV.flags.bAutoInc    = fread(fid,1,'ubit1');   % +bit c        increment filename after write
headers.ACV.flags.bAutoSave   = fread(fid,1,'ubit1');   % +bit d        autosave file after acquisition
headers.ACV.flags.bOverWrite  = fread(fid,1,'ubit1');   % +bit e        over write existing files

fread(fid,1,'ubit17'); % headers.ACV.flags.unused      =  bit 15-31    fill up to 32 bit - MUM 030627
headers.ACV.darate      = fread(fid,1,'uint32');        % 268 - 4 D/A clock frequency 
headers.ACV.dirate      = fread(fid,1,'uint32');        % 272 - 5 Pixel clock rate in mHz
headers.ACV.timeLeadX   = fread(fid,1,'int32');         % 276 - 6 DA Subsystems X time lead in microseconds 
% 280 - 7 pixels x and y 
headers.ACV.pixels.x    = fread(fid,1,'ubit16');        % was sxpix 
headers.ACV.pixels.y    = fread(fid,1,'ubit16');        % was sypix 
% 284 - 8
headers.ACV.clocksPerPixel.x  = fread(fid,1,'ubit16');  % was pppx 
headers.ACV.clocksPerPixel.y  = fread(fid,1,'ubit16');  % was pppy 
% 288 - 9 - scan parameters struct as replacement for old scan struct - MUM 030627
headers.ACV.scanParams.prePixels  = fread(fid,1,'ubit8');  % reclycled unused 'scan:16' - MUM 030715
headers.ACV.scanParams.postPixels = fread(fid,1,'ubit8');  % reclycled unused 'scan:16' - MUM 030715
headers.ACV.scanParams.yretrace   = fread(fid,1,'ubit16'); % not used any more with dtfu110.out and later - MUM 030627
% 292 - 10 + 11 
headers.ACV.num.images   = fread(fid,1,'ubit32');        % was 16 bit, extended to 32 - MUM 030627
headers.ACV.num.chans    = fread(fid,1,'ubit16');
headers.ACV.num.splitParts = fread(fid,1,'ubit16');      % Number of parts to split images into - MUM 030627
headers.ACV.timeLeadY    = fread(fid,1,'uint32');        % 300 - 12 - DA subsystems Y time lead in microseconds - MUM 030715
% 304 - 13
headers.ACV.shutterOffset.lines = fread(fid,1,'ubit16'); % Number of lines to open shutter in advance, minimum is 1 with dtfu121 and later
headers.ACV.shutterOffset.usec  = fread(fid,1,'ubit16'); % Same in microseconds for tighter control, for eventual use later
% 308 - 14 
headers.ACV.avg.num  = fread(fid,1,'ubit16');            % Number of images averaged/step 
fread(fid,1,'ubit16');            % headers.ACV.avg.toss = Number of images to toss before average - unused - MUM 030730
% 312 - 15 
headers.ACV.msize.x      = fread(fid,1,'ubit16');        % maximum x dimension, targetted was mxsize
headers.ACV.msize.y      = fread(fid,1,'ubit16');        % maximum y dimension, targetted was mysize
headers.ACV.scandur      = fread(fid,1,'single');        % 316 - 16 target scan duration (msec) 
headers.ACV.retracedur   = fread(fid,1,'single');        % 320 - 17 target retrace duration (msec) 
headers.ACV.ZStepSize    = fread(fid,1,'int32');         % 324 - 18 Z axis increment in nm 
headers.ACV.DwellTime    = fread(fid,1,'uint32');        % 328 - 19 Time in ms to wait during Z/T series scan 
fread(fid,1,'uint32');        % headers.ACV.cflags       = 332 - 20 config flag copies, now unused - MUM 031126
fread(fid,1,'uint32');        % headers.ACV.unused       = 336 - 21 - avg has been put to byte 14, because it is needed on the dt board - MUM 030721
headers.ACV.cycleTime    = fread(fid,1,'int32');         % 340 - 22 Time between acquisitions in TA mode 
headers.ACV.TacqCounter  = fread(fid,1,'int32');         % 344 - 23 number of cycles in TA mode 
headers.ACV.XStepSize    = fread(fid,1,'int32');         % 348 - 24 X axis increment in nm if XYZ (e.g. Sutter) stage used 
headers.ACV.YStepSize    = fread(fid,1,'int32');         % 352 - 25 Y axis increment in nm if XYZ (e.g. Sutter) stage used 
headers.ACV.vorlaufzeit  = fread(fid,1,'single');        % 356 - 26 target time for mirror movement in forward direction in ms
fread(fid,37,'int32');        % headers.ACV.rsvd         = 360 - 26 buffer to 38 longs (4 byte * 38 = 152 bytes) 
fread(fid,32,'ubit1');
% Third part - the ASV structure handles scan information
% fprintf('ASV: %i\n',ftell(fid))

% 0 Ranges for scan
headers.ASV.range.x     = fread(fid,1,'ubit16');
headers.ASV.range.y     = fread(fid,1,'ubit16');
% 4 Offsets for scan 
headers.ASV.offset.x    = fread(fid,1,'ubit16');
headers.ASV.offset.y    = fread(fid,1,'ubit16');
% 8 Parking position and orientation 
headers.ASV.scan.park   = fread(fid,1,'ubit16');
headers.ASV.scan.orient = fread(fid,1,'ubit16');
% 12 Flags 
fread(fid,1,'ubit32');         % headers.ASV.flags.rsvd  = all flags unused - MUM 030724
headers.ASV.zoom_fac    = fread(fid,1,'single');         % 16 Zoom factor 
fread(fid,9,'uint32');         % headers.ASV.reserv      = 20 reserved 
headers.ASV.Acm0.Gain   = fread(fid,1,'ubit16');         % 56 Analog module 0 gain, XPG
headers.ASV.Acm0.Offset = fread(fid,1,'ubit16');         % Analog module 0 offset, XPG 
headers.ASV.Acm1.Gain   = fread(fid,1,'ubit16');         % 60 Analog module 1 gain, XPG 
headers.ASV.Acm1.Offset = fread(fid,1,'ubit16');         % Analog module 1 offset, XPG 
fread(fid,48,'uint32');         % headers.ASV.rsvd       = 64 


% 4th part - WFM for waveform data is optional - if present, the structure is longer than 768 bytes
% fprintf('WFM: %i\n',ftell(fid))

%headers.WFM = fread(fid,48,'uint32');

fclose(fid);