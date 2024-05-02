function adaq2read(file)
%ADAQ2READ convert .adaq analog files to matlab format.
%
%   type adaq2read and knock yourself out ...
if ~nargin
    [ff,pth] = uigetfile('*.adaq','Choose ADAQ files to convert','multiselect','on');
    if ~iscell(ff)
        ff = {ff};
    end
    for n = 1 : length(ff)
        adaq2read([pth ff{n}])   
    end
    return
end

sep = strfind(file,'_');
outname = strrep(file(sep(end)+1:end),'.','_');
fid = fopen(file);
%  1   long   adaq_lFileSize;
%  2   long   adaq_lCheckSum_Header;            // beginning with the ilLenHeader field
%  3   long   adaq_lCheckSum_Data;
% 
%     // Common to all data files - complies with VDAQ Block file layout and values
%  4   long   adaq_lLenHeader;
%  5   long   adaq_lVersionID;
%  6   long   adaq_lFileType;                   // e.g. SUM, DC, RAW, ADAQ ...
%  7   long   adaq_lFileSubtype;                // e.g. FROM_VDAQ, FROM_ORA, FROM_ADAQ
%  8   long   adaq_lDataType;                   // e.g. DAT_... (uchar, ushort, long, float)
%  9   long   adaq_lSizeOf;                     // e.g. sizeof(short), sizeof(long)
% 
%     // information about the contents - ADAQ-Specific from here forward
%  10 long   adaq_lNChannels;                  // How many analog channels are collected
%  11 long   adaq_lNSamples;                   // How many samples per channel
%  12 long   adaq_lPreTriggerSamples;          // How many samples before the GO trigger
%  13 long   adaq_lNTrials;                    // Number of Trials for each condition (stored separately)
%  14 long   adaq_lNConditions;                // Number of Conditions (number of different stimuli)
%  15 long   adaq_lMaxSampleValue;             // 2047 for a 12-bit board
%  16 long   adaq_lSampleRateInHz;             // per channel sampling rate
x = fread(fid,16,'long');
sampling_rate = x(16);
pretrig = x(12);
nTrials = x(13);
nSamples = x(11);
nChannels = x(10);
nConditions = x(14);
headrLen = x(4);
%     float  adaq_lChMaxVoltage[ADAQ_MaxChannels = 16]; // Usually 10V, 5V, 0.5V or 0.05V
MaxVoltage = fread(fid,16,'float');
%     uchar  adaq_ucBipolar[ADAQ_MaxChannels = 16]; // 1 for bipolar (-max..max) input, 0 for unipolar (0..max)
Bipolar = fread(fid,16,'uchar');
%     uchar  adaq_ucListOfStimuli[256];        // The actual stimulus IDs
ListOfStimuli = fread(fid,256,'uchar');
%     char   adaq_sChComments[ADAQ_MaxChannels = 16][ADAQ_CommentSize = 16]; // Description of each channel's usage
ChComments = reshape(fread(fid,256,'char'),[16 16]);
%     uchar  adaq_ucTrialOrder[256];           // The Stim Indices, in order of presentation, for first 256 trials
TrialOrder = fread(fid,256,'uchar');
%     char   adaq_acRecordingDate[32];         // ctime() 26 chars
RecordingDate = fread(fid,256,'uchar');

% // Locate data
%     long   adaq_lCondOffs;  // location of Condition 0 (same as lLenHeader if adaq block is at start of file)
CondOstart = fread(fid,1,'long');
%     long   adaq_lChanSize;  // size of one channel's data:           lNSamples * lSizeof
ChanSize = fread(fid,1,'long');
%     long   adaq_lTrialSize; // size of one trial (all the channels): lNChannels * lChanSize
TrialSize = fread(fid,1,'long');
%     long   adaq_lCondSize;  // size of each Condition:               lNTrials * lTrialSize
condSize = fread(fid,1,'long');
%     // To read data for stim s, trial t, chan c, seek to
%     //    lSimOffs + s*lStimSize + t*lTrialSize + c*lChanSize
%     // and then read  lTrialSize  bytes.
%     // spare
%     char   _spare[128];
fclose(fid);
fid = fopen(file);
fread(fid,headrLen);
vars0 = '''sampling_rate'',''pretrig'',''nSamples'',''nTrials'',''nChannels'',''nConditions''';
vars = [];
for cond = 1 : nConditions
    for chan = 1 : nChannels
        x = fread(fid,nSamples,'short');
        eval(['Cond' num2str(cond-1) 'Trial0Chan' num2str(chan-1) ' = x;']); 
        vars = [vars ',' ['''Cond' num2str(cond-1) 'Trial0Chan' num2str(chan-1) '''']];
    end
end
eval(['save(''' outname ''',' vars0 vars ')'])

