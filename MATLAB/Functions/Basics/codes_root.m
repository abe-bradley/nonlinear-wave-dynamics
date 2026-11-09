function root_folder = codes_root()
%CODES_ROOT Repository folder containing the directory initializers.
%   Run init_dir first to make this helper available on the MATLAB path.
%   The root can have any name and is independent of the current directory.
%   Searches this file's ancestors for init_dir.m, so the helper can move
%   anywhere within the repository without changing its implementation.
    root_folder = fileparts(mfilename('fullpath'));
    while ~isfile(fullfile(root_folder,'init_dir.m'))
        parent_folder = fileparts(root_folder);
        if strcmp(parent_folder,root_folder) || isempty(parent_folder)
            error('codes_root:RootNotFound', ...
                'Cannot locate the project root: no ancestor folder contains init_dir.m.');
        end
        root_folder = parent_folder;
    end
end
