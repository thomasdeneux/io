function data = micam_read(fname)

% Input
if nargin<1
    fname = fn_getfile('*.rsh');
end

% Read file
data = get_micam_data(fname);


end


% class MICAMData(tuple):
%     """
%     Simple tuple wrapper, retrocompatible with previous
%     outputs of this module, but also brings name item getter
%     for better clarity of what data type is used in the code 
%     that use these outputs.
%     """
%     
%     def __new__(self,value):
%         return tuple.__new__(MICAMData,value)
%     
%     def keys(self):
%         return ["matrix","analog","diff_matrix"]
%     
%     def __getitem__(self,index):
%         if isinstance(index,str):     ###if index is a string, True
%             if index == "matrix":
%                 return super().__getitem__(0)
%             elif index == "analog":
%                 return super().__getitem__(1)
%             elif index == "diff_matrix":
%                 return super().__getitem__(2)
%             else :
%                 raise KeyError(f"Unknown key for MICAMData : {index} (valid keys are either 'matrix','analog','diff_matrix')")
%         else :
%             return super().__getitem__(index)
            

% def get_micam_data(InputPath):
% 
%     Files = read_rsh_filelist(InputPath)
%     Data, Signals, Rawdata= assemble_rsd_files(Files)
% 
%     return MICAMData((Data, Signals, Rawdata))

function micamData = get_micam_data(inputPath)
    files = read_rsh_filelist(inputPath);
    [data, signals, rawData] = assemble_rsd_files(files);
    micamData.matrix = data;
    micamData.analog = signals;
    micamData.diff_matrix = rawData;
end

    
% def read_rsh_filelist(InputPath):
% 
%     Files = []
%     with open(InputPath,'r') as F:
%         Line = F.readline()
%         while Line :
%             if "Data-File-List" in Line:
%                 FileLine = F.readline()
%                 while FileLine:
%                     if '.rsd' in FileLine or '.rsh' in FileLine :
%                         Files.append(os.path.join( os.path.dirname(InputPath) , FileLine.rstrip()))
%                     FileLine = F.readline()
%                 break
%             Line = F.readline()
%     #warnings.warn(f'Using RSD files : {Files}', category = RuntimeWarning ,stacklevel = 00)
%     #print(Files)
% 
%     return Files

function files = read_rsh_filelist(inputPath)
    files = {};
    fid = fopen(inputPath, 'r');
    line = fgetl(fid);
    while ischar(line)
        if contains(line, 'Data-File-List')
            fileLine = fgetl(fid);
            while ischar(fileLine)
                if contains(fileLine, {'.rsd', '.rsh'})
                    files{end+1} = fullfile(fileparts(inputPath), strtrim(fileLine));
                end
                fileLine = fgetl(fid);
            end
            break;
        end
        line = fgetl(fid);
    end
    fclose(fid);
end


% def read_rsd_data(InputPath):
%     
%     ### takes in a rsd file in argument, will transorm it in a (256, 100, 128) matrix
%     ### usable in python (numpy array)
% 
%     FramesPerFile = 256 #default behavior of RSD format
% 
%     msg = "Reading : " + os.path.split(InputPath)[1]
%     print(msg,end = "\r")
%     with open(InputPath ,'rb') as F:     ### rb stands for reading mode in binary mode
%         byte_content = F.read()
%     fmt = "<%dh" % (len(byte_content) // 2)
%     IntData = struct.unpack(fmt, byte_content)
%     _data = np.reshape(IntData,(FramesPerFile,100,128))
%     print( len(msg) * " ",end = "\r")
%     return _data

function data = read_rsd_data(inputPath)
    framesPerFile = 256;
    fid = fopen(inputPath, 'r');
    byteContent = fread(fid, 'int16');
    fclose(fid);
    data = reshape(byteContent, [128, 100, framesPerFile]);
    data = permute(data, [3, 2, 1]);
end


% def read_rsh_meta_data(InputPath):
%     import io, configparser
%     ini_str =  open(InputPath,"r").read() 
%     ini_str = "[section]\n acquisition_date"+ini_str.split("plsfile")[0].split("acquisition_date")[1]
%     ini_str = io.StringIO(ini_str)
%     config = configparser.ConfigParser()
%     config.read_file(ini_str)
%     return dict(config["section"])
% 
% def reshape_rsd_data(array3D, **kwargs):
%     Output = []
%     temp  = []
%     for I in range(np.shape(array3D)[0]):
%         for Line in range(int(np.shape(array3D)[1]/4)):
%             for J in range(np.shape(array3D)[2]):
%                 for K in range(4):
%                     temp.append(array3D[I,K+(Line*4),J])
%             Output.extend(temp)
%             temp  = []
%     if kwargs.get('reverse',False):
%         return np.asarray(Output)
%     else :
%         return np.invert(np.asarray(Output))

function reshapedData = reshape_rsd_data(array3D)
    [nx, ny4, nz] = size(array3D);
    ny = floor(ny4 / 4);
    reshapedData = array3D(:, 1:4*ny, :);
    reshapedData = reshape(...
        permute(...
            reshape(reshapedData, [nx, ny, 4, nz]), ...
            [1 3 2 4]), ...
        [nx ny*4, nz]);
end


% def assemble_rsd_files(InputPath):
%     

function [images, signals, rawData] = assemble_rsd_files(inputPaths)
    % #########
    % #configuration of the data inherent to the design of the RSD file as a "2D array".
    % #to store both the digital and analog signals as well as the images.
    % #it should not be necessary to ever change it.
    % 
    % ###to understand the next lines, micam ultima data format is necessary
    % #InputPath is the list of the paths to each the rsd file (16)

    % imagespan = [20,120]
    % SignalsSpan = [0,80]
    % FrameSpan = [10,12]
    % AIn1Span = [12,14]
    % AIn2Span = [14,16]
    % Stim1Span = [8,10]
    % Stim2Span = [6,8]
    % FilesPerSequence = 256

    framesPerFile = 256;
    imagespan = [21, 120];
    signalsSpan = [1, 80];
    frameSpan = [11, 12];
    aIn1Span = [13, 14];
    aIn2Span = [15, 16];
    stim1Span = [9, 10];
    stim2Span = [7, 8];
    
    % for ItemIndex in range(len(InputPath)):
    % 
    for itemIndex = 1:length(inputPaths)

        % Data = read_rsd_data(InputPath[ItemIndex])
        data = read_rsd_data(inputPaths{itemIndex});
        
        % if ItemIndex == 0 :
        %
        %     RawData = Data
        % 
        %     Image0 = Data[0,:,imagespan[0]:imagespan[1]]
        %     VarimagesImages = Data[1:,:,imagespan[0]:imagespan[1]]
        % 
        %     AnalogIn1 = Data[:,SignalsSpan[0]:SignalsSpan[1],AIn1Span[0]:AIn1Span[1]]
        %     AnalogIn2 = Data[:,SignalsSpan[0]:SignalsSpan[1],AIn2Span[0]:AIn2Span[1]]
        % 
        %     Stim1 = Data[:,SignalsSpan[0]:SignalsSpan[1],Stim1Span[0]:Stim1Span[1]]
        %     Stim2 = Data[:,SignalsSpan[0]:SignalsSpan[1],Stim2Span[0]:Stim2Span[1]]
        % 
        %     Frames = Data[:,SignalsSpan[0]:SignalsSpan[1],FrameSpan[0]:FrameSpan[1]]
        % 
        %     Images = np.empty((len(InputPath) * FilesPerSequence ,np.shape(VarimagesImages)[1] ,np.shape(VarimagesImages)[2]))
        % 
        %     for I in range(FilesPerSequence):
        % 
        %         if I == 0 :
        %             Images[I,:,:] = Image0
        %         else :
        %             Images[I,:,:] = VarimagesImages[I-1,:,:] + Image0      ###why + Image0 ??
        % 
        % else:
        % 
        %     RawData = np.append(RawData, Data, axis = 0)
        % 
        %     VarimagesImages = Data[:,:,imagespan[0]:imagespan[1]]
        % 
        %     AnalogIn1 = np.append(AnalogIn1 , Data[:,SignalsSpan[0]:SignalsSpan[1],AIn1Span[0]:AIn1Span[1]] , axis = 0)
        %     AnalogIn2 = np.append(AnalogIn2 , Data[:,SignalsSpan[0]:SignalsSpan[1],AIn2Span[0]:AIn2Span[1]] , axis = 0)
        % 
        %     Stim1 = np.append(Stim1 , Data[:,SignalsSpan[0]:SignalsSpan[1],Stim1Span[0]:Stim1Span[1]] , axis = 0)
        %     Stim2 = np.append(Stim2 , Data[:,SignalsSpan[0]:SignalsSpan[1],Stim2Span[0]:Stim2Span[1]] , axis = 0)
        % 
        %     Frames = np.append(Frames, Data[:,SignalsSpan[0]:SignalsSpan[1],FrameSpan[0]:FrameSpan[1]], axis = 0)
        % 
        %     for I in range(FilesPerSequence):
        % 
        %         Images[I+(ItemIndex * FilesPerSequence),:,:] = VarimagesImages[I,:,:] + Image0
        % 
        
        if itemIndex == 1
            rawData = data;
            image0 = data(1, :, imagespan(1):imagespan(2));
            varImages = data(2:end, :, imagespan(1):imagespan(2));
            analogIn1 = data(:, signalsSpan(1):signalsSpan(2), aIn1Span(1):aIn1Span(2));
            analogIn2 = data(:, signalsSpan(1):signalsSpan(2), aIn2Span(1):aIn2Span(2));
            stim1 = data(:, signalsSpan(1):signalsSpan(2), stim1Span(1):stim1Span(2));
            stim2 = data(:, signalsSpan(1):signalsSpan(2), stim2Span(1):stim2Span(2));
            frames = data(:, signalsSpan(1):signalsSpan(2), frameSpan(1):frameSpan(2));
            images = zeros([length(inputPaths) * framesPerFile, size(varImages, 2), size(varImages, 3)]);
            images(1, :, :) = image0;
            images(2:framesPerFile, :, :) = varImages + image0;
        else
            rawData = cat(1, rawData, data);
            varImages = data(:, :, imagespan(1):imagespan(2));
            analogIn1 = cat(1, analogIn1, data(:, signalsSpan(1):signalsSpan(2), aIn1Span(1):aIn1Span(2)));
            analogIn2 = cat(1, analogIn2, data(:, signalsSpan(1):signalsSpan(2), aIn2Span(1):aIn2Span(2)));
            stim1 = cat(1, stim1, data(:, signalsSpan(1):signalsSpan(2), stim1Span(1):stim1Span(2)));
            stim2 = cat(1, stim2, data(:, signalsSpan(1):signalsSpan(2), stim2Span(1):stim2Span(2)));
            frames = cat(1, frames, data(:, signalsSpan(1):signalsSpan(2), frameSpan(1):frameSpan(2)));
            images((itemIndex-1)*framesPerFile + (1:framesPerFile), :, :) = varImages + image0;
        end
    end
    
    % print("Processing.",end = "\r")
    % AnalogIn1 = rescale_micam_signal(reshape_rsd_data(AnalogIn1))
    % AnalogIn2 = rescale_micam_signal(reshape_rsd_data(AnalogIn2))
    % print("Processing..",end = "\r")
    % Stim1 = rescale_micam_signal(reshape_rsd_data(Stim1))
    % Stim2 = rescale_micam_signal(reshape_rsd_data(Stim2))
    % 
    % Frames = reshape_rsd_data(Frames,reverse = True)
    % print("Processing...",end = "\r")
    % RawData = np.moveaxis(RawData,0,2).astype(np.int16)
    % Images = np.moveaxis(Images,0,2).astype(np.int16)
    % 
    % Signals = {'AI1' : AnalogIn1, 'AI2' : AnalogIn2, 'Stim1' : Stim1, 'Stim2' : Stim2, 'FrameAcq' : Frames}
    % print(" "*16,end = "\r")
    % return Images, Signals, RawData

    analogIn1 = rescale_micam_signal(reshape_rsd_data(analogIn1));
    analogIn2 = rescale_micam_signal(reshape_rsd_data(analogIn2));
    stim1 = rescale_micam_signal(reshape_rsd_data(stim1));
    stim2 = rescale_micam_signal(reshape_rsd_data(stim2));
    frames = reshape_rsd_data(frames);
    rawData = permute(int16(rawData), [2, 3, 1]);
    images = permute(int16(images), [2, 3, 1]);
    
    signals.AI1 = analogIn1;
    signals.AI2 = analogIn2;
    signals.Stim1 = stim1;
    signals.Stim2 = stim2;
    signals.FrameAcq = frames;
end


% def rescale_micam_signal(signalist, upper_voltage = 5):
%     #signal between 0 and 5 V (or upper_voltage if provided)
%     return ((signalist.astype('float32') * (upper_voltage/56500)) + ((upper_voltage/56500) * 23730) ).astype('float16')

function rescaledSignal = rescale_micam_signal(signalList)
    upperVoltage = 5;
    rescaledSignal = ((double(signalList) * (upperVoltage/56500)) + ((upperVoltage/56500) * 23730));
end



% def get_micam_trigger_offset(signalist):
%     sys.path.append(os.path.dirname(os.path.dirname(__file__)))
%     from sigprocess import measurements
%     #signalist = rescale_micam_signal(signalist)
%     dictsignal = measurements.DetectHighPulses(signalist, 3, 1, 80000)
%     return 500 - ( dictsignal['count'] - 4 )

% if __name__ == "__main__" :
%     root = r"D:\DATA\MICAM\VSD\Mouse63\210521_VSD1"
% 
%     inpath = "Behavior-1-1.rsh"
%     DATA = get_micam_data(os.path.join(root,inpath))
%     
%     real_vsd_data = DATA[0]
%     #or
%     real_vsd_data= DATA["matrix"]
%     
%     analog_signals = DATA[1]
%     #or
%     analog_signals = DATA["analog"]
%     print(analog_signals.keys())
%     #>> 'AI1' : AnalogIn1, 'AI2' : AnalogIn2, 'Stim1' : Stim1, 'Stim2' : Stim2, 'FrameAcq' : Frames incremental timings
%     
%     
%     delta0_vsd_data = DATA[2]
%     #or
%     delta0_vsd_data = DATA["diff_matrix"]
    
    
    
