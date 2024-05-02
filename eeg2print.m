function eeg2print(exprs,k)

% this function provide input to adq2vu to print all the elecrical data
% automaticall
%INPUTS: exprs:experiments. e.g. [1 5 7] or [2 3] or 2:5
%        %k  the number of conditions 
exprs = exprs(:)';
for conds = 1 : k
    for Es = exprs
        files = struct2cell(dir(['E' num2str(Es) '*_adaq.mat']));
        files = files(1,:);
        [~,ord] = sort(cellfun(@(x) str2num(strrep(x(4:end),'_adaq.mat','')),files));
        files = files(ord);
        adq2vu(1,conds-1,[4 200],[1 2500],0,1,files)
        close
    end
end