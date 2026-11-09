function root_folder = codes_root()
%CODES_ROOT Folder containing this helper and the directory initializers.
%   Run init_dir first to make this helper available on the MATLAB path.
%   The root can have any name and is independent of the current directory.
    root_folder = fileparts(mfilename('fullpath'));
end
