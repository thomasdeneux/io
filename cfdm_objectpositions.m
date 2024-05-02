function varargout = cfdm_objectpositions(f)
% function [FIGPOS POS] = cfdm_objectpositions
% function pos = cfdm_objectpositions(f)
% function [FIGPOS POS] = cfdm_objectpositions('frame')
% function cfdm_objectpositions('design')

USEMAT = 1;

if USEMAT
    
    load([fn_cd('2ps') '/save/cfdm_objectpositions'])
    
else
    
    FIGPOS = [40    40  1192   688];
    POS.ha1 = [24  416  256  256];
    POS.ha2 = [328  416  256  256];
    POS.hu1 = [56  384  192   16];
    POS.hu2 = [360  384  192   16];
    POS.slideraxes = [600  416    8  256];
    POS.htime = [40    24  1088   304];
    POS.listdir = [836  408  164  268];
    POS.listfile = [1008   408   180   268];
    
end

if nargin==1 && strcmp(f,'frame')
    
    varargout = {FIGPOS POS};
    return
    
elseif (nargin==0 && nargout==0) || (nargin==1 && strcmp(f,'design'))

    cmd = fn_framedesign(FIGPOS,POS);
    clear FIGPOS POS, eval(cmd')
    save([fn_cd('2ps') '/save/cfdm_objectpositions'],'FIGPOS','POS','cmd')
    varargout = {};
    return

end

axname = {'ha1','ha2','htime','slideraxes'};
axcorr = FIGPOS([3 4 3 4]);
for k=1:length(axname)
    g = axname{k};
    POS.(g) = POS.(g) ./ axcorr;
end

if nargin==0
    varargout = {FIGPOS POS};
else
    varargout = {POS.(f)};
end

