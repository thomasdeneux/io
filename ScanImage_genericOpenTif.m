function [Aout,header] = ScanImage_genericOpenTif(filename,varargin)
% GENERICOPENTIF   - Opens an Image file with TIF extension.
%   GENERICOPENTIF opens a TIF file and store its contents as array Aout.
%   Filename is the file name with extension. 
%   
%   VARARGIN are paramameter value pairs that specify the output type.
%   Possible values include:
%
%    'filter'                1 or 0      Apply blocksize x blocksize Median Filter
%    'blocksize'             > 1         Blocksize for filter.
%    'splitIntoCellArray'    1 or 0      Output each input channel in cell array.
%    'linescan'              1 or 0      Reshape output into single frame by concatentating
%                                             to the bottom of the image.
%
%   See also CONVERTSTACKTOLS, PARSEHEADER

Aout=[];
header=[];

% Parse the inputs....
filter=0;
blocksize=3;
splitIntoCellArray=0;
linescan=0;
if nargin > 1
    % Parse input parameter pairs and rewrite values.
    counter=1;
    while counter+1 <= length(varargin)
        eval([varargin{counter} '=[(varargin{counter+1})];']);
        counter=counter+2;
    end
end


h = waitbar(0,'Opening Tif image...', 'Name', 'Open TIF Image', 'Pointer', 'watch');
try
    info=imfinfo(filename);
    frames = length(info);
    header=info(1).ImageDescription;
    header=parseHeader(header);
    for i = 1:frames
        waitbar(i/frames,h, ['Loading Frame Number ' num2str(i)]);    
        Aout(:,:,i) = imread(filename, i);
		if filter
			Aout(:,:,i)=medfilt2(Aout(:,:,i),[blocksize blocksize]);
		end
    end
    waitbar(1,h, 'Done');
    close(h);
catch
    close(h);
    disp(['Cant load file: ' filename ]);
end


% Pushes the data into cell arrays according to the number of channels read....
if splitIntoCellArray
    channels=header.acq.numberOfChannelsAcquire;
    for channelCounter=1:channels
        data{channelCounter}=Aout(:,:,channelCounter:channels:end);
    end
    Aout=data;
end

if linescan
    if iscell(Aout)
        for j=1:length(Aout)
            Aout{j}=convertStackToLS(Aout{j});
        end
    else
        Aout=convertStackToLS(Aout);
    end
end
