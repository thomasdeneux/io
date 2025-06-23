function data = rsd_read(fname)

% Input
if nargin<1
    fname = fn_getfile('*.rsd');
end

% Read file
data = get_micam_data(fname);


end


function micamData = get_micam_data(inputPath)
    files = read_rsh_filelist(inputPath);
    [data, signals, rawData] = assemble_rsd_files(files);
    micamData.matrix = data;
    micamData.analog = signals;
    micamData.diff_matrix = rawData;
end

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

function data = read_rsd_data(inputPath)
    framesPerFile = 256;
    fid = fopen(inputPath, 'r');
    byteContent = fread(fid, 'int16');
    fclose(fid);
    data = reshape(byteContent, [100, 128, framesPerFile]);
    data = permute(data, [3, 1, 2]);
end

function reshapedData = reshape_rsd_data(array3D)
    output = [];
    for i = 1:size(array3D, 1)
        for line = 1:floor(size(array3D, 2) / 4)
            for j = 1:size(array3D, 3)
                temp = array3D(i, (line-1)*4 + (1:4), j);
                output = [output; temp(:)'];
            end
        end
    end
    reshapedData = bitcmp(output, 'uint16');
end

function [images, signals, rawData] = assemble_rsd_files(inputPaths)
    framesPerFile = 256;
    imagespan = [21, 121];
    signalsSpan = [1, 81];
    frameSpan = [11, 13];
    aIn1Span = [13, 15];
    aIn2Span = [15, 17];
    stim1Span = [9, 11];
    stim2Span = [7, 9];
    
    for itemIndex = 1:length(inputPaths)
        data = read_rsd_data(inputPaths{itemIndex});
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

function rescaledSignal = rescale_micam_signal(signalList)
    upperVoltage = 5;
    rescaledSignal = ((double(signalList) * (upperVoltage/56500)) + ((upperVoltage/56500) * 23730));
end
