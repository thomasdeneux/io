# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:47:27 2020

@author: Timothe
"""

import struct
import numpy as np
import os, sys




class MICAMData(tuple):
    """
    Simple tuple wrapper, retrocompatible with previous
    outputs of this module, but also brings name item getter
    for better clarity of what data type is used in the code 
    that use these outputs.
    """
    
    def __new__(self,value):
        return tuple.__new__(MICAMData,value)
    
    def keys(self):
        return ["matrix","analog","diff_matrix"]
    
    def __getitem__(self,index):
        if isinstance(index,str):     ###if index is a string, True
            if index == "matrix":
                return super().__getitem__(0)
            elif index == "analog":
                return super().__getitem__(1)
            elif index == "diff_matrix":
                return super().__getitem__(2)
            else :
                raise KeyError(f"Unknown key for MICAMData : {index} (valid keys are either 'matrix','analog','diff_matrix')")
        else :
            return super().__getitem__(index)
        
# def get_rsh_files(SessionPath):
    
#     Files = strings.RegFileSearch(SessionPath,'.*\.rsh$')
#     Files = strings.AlphaNum_Sort(Files)

#     return Files



        
def get_micam_data(InputPath):

    Files = read_rsh_filelist(InputPath)
    Data, Signals, Rawdata= assemble_rsd_files(Files)

    return MICAMData((Data, Signals, Rawdata))

def read_rsh_filelist(InputPath):

    Files = []
    with open(InputPath,'r') as F:
        Line = F.readline()
        while Line :
            if "Data-File-List" in Line:
                FileLine = F.readline()
                while FileLine:
                    if '.rsd' in FileLine or '.rsh' in FileLine :
                        Files.append(os.path.join( os.path.dirname(InputPath) , FileLine.rstrip()))
                    FileLine = F.readline()
                break
            Line = F.readline()
    #warnings.warn(f'Using RSD files : {Files}', category = RuntimeWarning ,stacklevel = 00)
    #print(Files)

    return Files

def read_rsd_data(InputPath):
    
    ### takes in a rsd file in argument, will transorm it in a (256, 100, 128) matrix
    ### usable in python (numpy array)

    FramesPerFile = 256 #default behavior of RSD format

    msg = "Reading : " + os.path.split(InputPath)[1]
    print(msg,end = "\r")
    with open(InputPath ,'rb') as F:     ### rb stands for reading mode in binary mode
        byte_content = F.read()
    fmt = "<%dh" % (len(byte_content) // 2)
    IntData = struct.unpack(fmt, byte_content)
    _data = np.reshape(IntData,(FramesPerFile,100,128))
    print( len(msg) * " ",end = "\r")
    return _data

def read_rsh_meta_data(InputPath):
    import io, configparser
    ini_str =  open(InputPath,"r").read() 
    ini_str = "[section]\n acquisition_date"+ini_str.split("plsfile")[0].split("acquisition_date")[1]
    ini_str = io.StringIO(ini_str)
    config = configparser.ConfigParser()
    config.read_file(ini_str)
    return dict(config["section"])

def reshape_rsd_data(array3D, **kwargs):
    Output = []
    temp  = []
    for I in range(np.shape(array3D)[0]):
        for Line in range(int(np.shape(array3D)[1]/4)):
            for J in range(np.shape(array3D)[2]):
                for K in range(4):
                    temp.append(array3D[I,K+(Line*4),J])
            Output.extend(temp)
            temp  = []
    if kwargs.get('reverse',False):
        return np.asarray(Output)
    else :
        return np.invert(np.asarray(Output))

def assemble_rsd_files(InputPath):
    #########
    #configuration of the data inherent to the design of the RSD file as a "2D array".
    #to store both the digital and analog signals as well as the images.
    #it should not be necessary to ever change it.
    
    ###to understand the next lines, micam ultima data format is necessary
    #InputPath is the list of the paths to each the rsd file (16)
    
    
    imagespan = [20,120]
    SignalsSpan = [0,80]
    FrameSpan = [10,12]
    AIn1Span = [12,14]
    AIn2Span = [14,16]
    Stim1Span = [8,10]
    Stim2Span = [6,8]
    FilesPerSequence = 256
    ###############
    
    for ItemIndex in range(len(InputPath)):

        Data = read_rsd_data(InputPath[ItemIndex])
        
        if ItemIndex == 0 :

            RawData = Data

            Image0 = Data[0,:,imagespan[0]:imagespan[1]]
            VarimagesImages = Data[1:,:,imagespan[0]:imagespan[1]]

            AnalogIn1 = Data[:,SignalsSpan[0]:SignalsSpan[1],AIn1Span[0]:AIn1Span[1]]
            AnalogIn2 = Data[:,SignalsSpan[0]:SignalsSpan[1],AIn2Span[0]:AIn2Span[1]]

            Stim1 = Data[:,SignalsSpan[0]:SignalsSpan[1],Stim1Span[0]:Stim1Span[1]]
            Stim2 = Data[:,SignalsSpan[0]:SignalsSpan[1],Stim2Span[0]:Stim2Span[1]]

            Frames = Data[:,SignalsSpan[0]:SignalsSpan[1],FrameSpan[0]:FrameSpan[1]]

            Images = np.empty((len(InputPath) * FilesPerSequence ,np.shape(VarimagesImages)[1] ,np.shape(VarimagesImages)[2]))

            for I in range(FilesPerSequence):

                if I == 0 :
                    Images[I,:,:] = Image0
                else :
                    Images[I,:,:] = VarimagesImages[I-1,:,:] + Image0      ###why + Image0 ??

        else:

            RawData = np.append(RawData, Data, axis = 0)

            VarimagesImages = Data[:,:,imagespan[0]:imagespan[1]]

            AnalogIn1 = np.append(AnalogIn1 , Data[:,SignalsSpan[0]:SignalsSpan[1],AIn1Span[0]:AIn1Span[1]] , axis = 0)
            AnalogIn2 = np.append(AnalogIn2 , Data[:,SignalsSpan[0]:SignalsSpan[1],AIn2Span[0]:AIn2Span[1]] , axis = 0)

            Stim1 = np.append(Stim1 , Data[:,SignalsSpan[0]:SignalsSpan[1],Stim1Span[0]:Stim1Span[1]] , axis = 0)
            Stim2 = np.append(Stim2 , Data[:,SignalsSpan[0]:SignalsSpan[1],Stim2Span[0]:Stim2Span[1]] , axis = 0)

            Frames = np.append(Frames, Data[:,SignalsSpan[0]:SignalsSpan[1],FrameSpan[0]:FrameSpan[1]], axis = 0)

            for I in range(FilesPerSequence):

                Images[I+(ItemIndex * FilesPerSequence),:,:] = VarimagesImages[I,:,:] + Image0
                
                
    print("Processing.",end = "\r")
    AnalogIn1 = rescale_micam_signal(reshape_rsd_data(AnalogIn1))
    AnalogIn2 = rescale_micam_signal(reshape_rsd_data(AnalogIn2))
    print("Processing..",end = "\r")
    Stim1 = rescale_micam_signal(reshape_rsd_data(Stim1))
    Stim2 = rescale_micam_signal(reshape_rsd_data(Stim2))
    
    Frames = reshape_rsd_data(Frames,reverse = True)
    print("Processing...",end = "\r")
    RawData = np.moveaxis(RawData,0,2).astype(np.int16)
    Images = np.moveaxis(Images,0,2).astype(np.int16)

    Signals = {'AI1' : AnalogIn1, 'AI2' : AnalogIn2, 'Stim1' : Stim1, 'Stim2' : Stim2, 'FrameAcq' : Frames}
    print(" "*16,end = "\r")
    return Images, Signals, RawData

def rescale_micam_signal(signalist, upper_voltage = 5):
    #signal between 0 and 5 V (or upper_voltage if provided)
    return ((signalist.astype('float32') * (upper_voltage/56500)) + ((upper_voltage/56500) * 23730) ).astype('float16')

def get_micam_trigger_offset(signalist):
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from sigprocess import measurements
    #signalist = rescale_micam_signal(signalist)
    dictsignal = measurements.DetectHighPulses(signalist, 3, 1, 80000)
    return 500 - ( dictsignal['count'] - 4 )

if __name__ == "__main__" :
    root = r"E:\Users\Thomas\LocalData\2025-02 Noe\StimC2_S1"

    inpath = "NoStim_C2S1_-6-1.rsh"
    DATA = get_micam_data(os.path.join(root,inpath))
    
    real_vsd_data = DATA[0]
    #or
    real_vsd_data= DATA["matrix"]
    
    analog_signals = DATA[1]
    #or
    analog_signals = DATA["analog"]
    print(analog_signals.keys())
    #>> 'AI1' : AnalogIn1, 'AI2' : AnalogIn2, 'Stim1' : Stim1, 'Stim2' : Stim2, 'FrameAcq' : Frames incremental timings
    
    
    delta0_vsd_data = DATA[2]
    #or
    delta0_vsd_data = DATA["diff_matrix"]
    
    
    
    
    
