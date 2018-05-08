% This script organizes relevant files into a certain logical order and
% passes them onto the blob analyze function to calculate average value of
% each channel within the structures defined by the mask image.

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

if isempty(file_search('blb_\w+.txt',folder))
    mkdir(folder,'Average Images')
    col_labels = blob_analyze(imageset,sizemin,sizemax,exp_name,Achannel,folder);
end
addpath(fullfile(folder,'Average Images'))