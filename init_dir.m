clear
clc

%% Initialize this file's root, excluding packages, recycling, and scratch exports
codes_folder = fileparts(mfilename('fullpath'));

code_folders = strsplit(genpath(codes_folder),pathsep);
recycling_folder = fullfile(codes_folder,'Recycling');
recycling_subfolder = [recycling_folder filesep];
is_recycling_folder = strcmpi(code_folders,recycling_folder) | ...
    strncmpi(code_folders,recycling_subfolder,length(recycling_subfolder));

packages_folder = fullfile(codes_folder,'Packages');
packages_subfolder = [packages_folder filesep];
is_packages_folder = strcmpi(code_folders,packages_folder) | ...
    strncmpi(code_folders,packages_subfolder,length(packages_subfolder));

% Scratch exports can contain stale copies of live functions and libraries.
temporary_folder = fullfile(codes_folder,'.temp');
temporary_subfolder = [temporary_folder filesep];
is_temporary_folder = strcmpi(code_folders,temporary_folder) | ...
    strncmpi(code_folders,temporary_subfolder,length(temporary_subfolder));

old_paths = strsplit(path,pathsep);
is_code_path = strcmpi(old_paths,codes_folder) | ...
    startsWith(old_paths,[codes_folder filesep],'IgnoreCase',true);
if any(is_code_path)
    rmpath(old_paths{is_code_path});
end
cd(codes_folder);
addpath(strjoin(code_folders(~is_recycling_folder & ~is_packages_folder ...
    & ~is_temporary_folder),pathsep));
clearvars -except codes_folder

fprintf("\nDirectory initialized!\n\n")
