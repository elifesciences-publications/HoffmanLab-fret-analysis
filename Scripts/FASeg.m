% This script deals with only segmentation protocols for focal adhesion
% structures. It generates focal adhesion masks based on three input
% parameters or "blob_params"

rehash
if isempty(file_search('fa_\w+.TIF',folder))
    mkdir(folder,'FA Images')
    fa_gen(['bsa_' prefix exp_name '\w+' Achannel '.TIF'],blob_params,folder)
end
addpath(fullfile(folder,'FA Images'))