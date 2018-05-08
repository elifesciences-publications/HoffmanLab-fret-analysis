% This script applies masks generated on focal adhesions to all relevant
% imaging channels.
% For corrected FRET, this includes "bsa","bsd", and "cna" channels
% For FRET efficiency, it adds on "eff" and "dpa" channels
% For calibrated FRET efficiency, it adds on "for" and "ext" channels

rehash
imageset = {['bsa_' prefix exp_name '\w+' Achannel '.TIF'],...
    ['bsd_' prefix exp_name '\w+' Dchannel '.TIF']};
if strcmpi(FRETeff,'y')
    imageset{end+1} = ['dpa_' prefix exp_name '\w+' FRETchannel '.TIF'];
end
imageset{end+1} = ['cna_' prefix exp_name '\w+' FRETchannel '.TIF'];
if strcmpi(FRETeff,'y')
    imageset{end+1} = ['eff_' prefix exp_name '\w+' FRETchannel '.TIF'];
end
if strcmpi(Force,'y')
    imageset{end+1} = ['force_' prefix exp_name '\w+' FRETchannel '.TIF'];
    imageset{end+1} = ['ext_' prefix exp_name '\w+' FRETchannel '.TIF'];
end
imageset{end+1} = ['fa_bsa_' prefix exp_name '\w+' Achannel '.TIF'];


if isempty(file_search('masked_\w+.TIF',folder))
    mkdir(folder,'Masked Images')
    app_mask(imageset,Achannel,folder)
end
addpath(fullfile(folder,'Masked Images'))