% This script figures out what channel to draw cell boundaries on,
% and calls a function to manually draw cell masks.

rehash

if ~exist('Achannel','var')
    Achannel = '';
end

if strcmpi(BoundaryChannel,Achannel)
    files = ['bsa_' prefix exp_name '\w+' Achannel '.TIF'];
elseif strcmpi(BoundaryChannel,FRETchannel)
    files = ['c_' prefix exp_name '\w+' FRETchannel '.TIF'];
elseif strcmpi(BoundaryChannel,Dchannel)
    files = ['bsd_' prefix exp_name '\w+' Dchannel '.TIF'];
end

if isempty(file_search('polymask_\w+',folder))
    mkdir(fullfile(folder,'Cell Mask Images'))
    cell_outline_manual(files,folder);
end

addpath(fullfile(folder,'Cell Mask Images'))