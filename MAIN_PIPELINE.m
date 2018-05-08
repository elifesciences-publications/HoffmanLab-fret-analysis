function MAIN_PIPELINE(varargin)

%% Header
% This is a main function for running FRET analysis on a single
% experimental group. The Exp_Params text file in the folder containing the data
% should describe all relevant information for a particular experimental design.
% Sample Exp_params text files can be found in the GitHub repository.

%% Parameter Processing
% Check to see if anything has been passed as parameter, if anything has
% been passed, make sure it is a folder.
addpath(genpath(pwd))
if (not(isempty(varargin)))
    if (exist(varargin{1},'dir'))
        folder = varargin{1};
    else
        error('Expected first parameter to be a folder with images to process.');
    end
else
    %% Set up
    clear;
    close all;
    clc;
end

%% Read in Pre-processing parameters

if (exist('folder','var'))
    [~,params_file] = GetParamsFile(folder); %#ok<ASGLU>
else
    [folder,params_file] = GetParamsFile; %#ok<ASGLU>
end
ProcessParamsFile;

%% Pre-process
PreProcessImgs

%% Run FRET correction, does FRET efficiency if desired
FRET

%% Segmentation
if strcmpi(segmentation,'y') && strcmpi(structure,'FA')
    FASeg
end

%% Mask Images
if strcmpi(mask,'y') && strcmpi(structure,'FA') && strcmpi(segmentation,'y')
    FAMask
end

%% Blob Analysis
if strcmpi(segmentation,'y') && strcmpi(banalyze,'y') && strcmpi(structure,'FA')
    BlobAnalyze
end

%% Draw Boundaries
if strcmpi(draw_boundaries,'y')
    DrawBoundaries
end

%% Add Boundary Properties to Blob Analysis Results
if strcmpi(add_boundary_props,'y')
    AddBoundaryProps
end

%% Update parameters
SaveParams