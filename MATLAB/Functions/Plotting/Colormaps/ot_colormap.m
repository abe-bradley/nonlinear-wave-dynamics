function [colors,info] = ot_colormap(name,N,varargin)
%OT_COLORMAP Native MATLAB and attributed external palettes through one API.
%   C = OT_COLORMAP(NAME,N) returns an N-by-3 RGB array. N defaults to 256.
%   Names are case insensitive. Original slanCM numeric IDs are retained;
%   removed IDs are not reassigned. Native names use installed MATLAB maps.
%   [C,INFO] also returns the palette's source and license metadata.
%   C = OT_COLORMAP(NAME,N,'Sampling','colorbrewer') selects official
%   count-specific ColorBrewer colors. Other calls retain existing sampling.
%
%   Examples:
%       colormap(ax,ot_colormap('thermal-2',256))
%       line_colors = ot_colormap('guppy',10);
%
%   External RGB data is extracted from slanCM and retains its collection
%   and original-source notices. See COLORMAPS.txt, THIRD_PARTY_NOTICES.txt,
%   and LICENSES. This function does not require slanCM installed.

    narginchk(1,4)
    if nargin < 2
        N = 256;
    end
    if ~isnumeric(N) || ~isreal(N) || ~isscalar(N) || ~isfinite(N) ...
            || N < 0 || N ~= fix(N)
        error('ot_colormap:InvalidSize', ...
            'The palette length must be a finite, nonnegative integer.')
    end
    N = double(N);
    sampling = 'default';
    if ~isempty(varargin)
        if numel(varargin) ~= 2 || ~scalar_text(varargin{1}) ...
                || ~strcmpi(strtrim(char(varargin{1})),'Sampling') ...
                || ~scalar_text(varargin{2})
            error('ot_colormap:InvalidOptions', ...
                'Use the name-value option ''Sampling'',''colorbrewer'' or ''default''.')
        end
        sampling = lower(strtrim(char(varargin{2})));
        if ~ismember(sampling,{'default','colorbrewer'})
            error('ot_colormap:InvalidSampling','Unknown sampling policy "%s".',sampling)
        end
    end

    persistent catalog palettes
    resource_folder = fullfile(fileparts(mfilename('fullpath')),'data');
    if isempty(catalog)
        catalog = jsondecode(fileread(fullfile(resource_folder, ...
            'palette_catalog.json')));
    end
    entries = catalog.entries;
    if isnumeric(name)
        if ~isreal(name) || ~isscalar(name) || ~isfinite(name) ...
                || name < 1 || name ~= fix(name)
            error('ot_colormap:InvalidID', ...
                'A numeric palette ID must be a positive integer.')
        end
        index = find([entries.id] == double(name),1);
        if isempty(index)
            error('ot_colormap:UnknownPaletteID', ...
                'Palette ID %g is not in the published catalog.',double(name))
        end
    elseif (ischar(name) && isrow(name)) ...
            || (isstring(name) && isscalar(name) && ~ismissing(name))
        name = strtrim(char(name));
        index = find(strcmpi({entries.name},name),1);
        if isempty(index)
            error('ot_colormap:UnknownPalette', ...
                'Palette "%s" is not in the published catalog.',name)
        end
    else
        error('ot_colormap:InvalidName', ...
            'Specify one palette name or an original numeric palette ID.')
    end
    info = entries(index);
    info.sampling = sampling;
    if strcmp(sampling,'colorbrewer')
        [colors,info] = colorbrewer_colors(info,N,resource_folder);
        return
    end

    switch info.provider
        case 'matlab'
            % Explicit native registry; do not execute arbitrary function
            % names, or a project file shadowing a MATLAB palette function.
            native_names = {'parula','turbo','hsv','hot','cool','spring', ...
                'summer','autumn','winter','gray','bone','copper','pink', ...
                'jet','lines','colorcube','prism','flag','white', ...
                'abyss','sky','nebula'};
            if ~ismember(info.function_name,native_names)
                error('ot_colormap:InvalidCatalog','Invalid native palette entry.')
            end
            native_file = which(info.function_name);
            if isempty(native_file)
                error('ot_colormap:NativeUnavailable', ...
                    'Native palette "%s" is unavailable in MATLAB %s.', ...
                    info.name,version('-release'))
            end
            if ~startsWith(native_file,[matlabroot filesep],'IgnoreCase',ispc)
                error('ot_colormap:NativeShadowed', ...
                    'Palette "%s" resolves outside MATLAB: %s. Remove that path conflict.', ...
                    info.name,native_file)
            end
            colors = feval(info.function_name,N);
            info.matlab_release = version('-release');
        case 'external'
            if isempty(palettes)
                loaded = load(fullfile(resource_folder, ...
                    'external_palettes.mat'),'palettes');
                palettes = loaded.palettes;
            end
            if ~isfield(palettes,info.data_key)
                error('ot_colormap:MissingPaletteData', ...
                    'External data for "%s" is missing.',info.name)
            end
            base_colors = palettes.(info.data_key);
            if N == 0
                colors = zeros(0,3);
            else
                % Retain slanCM's linear interpolation of 256 base colors,
                % including its final base color when N is one.
                colors = interp1(1:size(base_colors,1),base_colors, ...
                    linspace(1,size(base_colors,1),N),'linear');
            end
            info.matlab_release = '';
        otherwise
            error('ot_colormap:InvalidCatalog','Unknown palette provider.')
    end
end

function valid = scalar_text(value)
    valid = (ischar(value) && isrow(value)) ...
        || (isstring(value) && isscalar(value) && ~ismissing(value));
end

function [colors,info] = colorbrewer_colors(info,N,resource_folder)
% Independent sampler using original ColorBrewer count-specific RGB data.
    persistent source
    if isempty(source)
        source = jsondecode(fileread(fullfile(resource_folder, ...
            'colorbrewer_schemes.json')));
    end
    index = find(strcmpi({source.schemes.name},info.name),1);
    if isempty(index)
        error('ot_colormap:NotColorBrewer', ...
            'ColorBrewer sampling is available only for its 35 scheme names; "%s" is not one.', ...
            info.name)
    end
    scheme = source.schemes(index);
    variants = scheme.variants;
    counts = [variants.count];
    exact = find(counts == N,1);
    [~,largest] = max(counts);
    base = double(variants(largest).colors)/255;
    if N == 0
        colors = zeros(0,3);
        method = 'empty';
    elseif ~isempty(exact)
        colors = double(variants(exact).colors)/255;
        method = 'designed';
    elseif strcmp(scheme.type,'qualitative')
        % Small requests use the leading category colors. Excess categories
        % repeat the largest designed set, without blending new colors.
        colors = base(1 + mod(0:N-1,size(base,1)),:);
        method = 'category-selection-or-repeat';
    elseif N < min(counts)
        % ColorBrewer defines no one/two-color selections. Use the center or
        % endpoints of its smallest designed selection as an explicit extension.
        [~,smallest] = min(counts);
        small = double(variants(smallest).colors)/255;
        if N == 1
            colors = small(ceil(size(small,1)/2),:);
        else
            colors = small([1 end],:);
        end
        method = 'small-count-extension';
    else
        knots = 1:size(base,1);
        query = linspace(knots(1),knots(end),N);
        if strcmp(scheme.type,'diverging')
            % Treat the two arms separately to retain the neutral center.
            middle = ceil(size(base,1)/2);
            left = query <= middle;
            colors = zeros(N,3);
            colors(left,:) = interp1(knots(1:middle),base(1:middle,:), ...
                query(left),'pchip');
            colors(~left,:) = interp1(knots(middle:end),base(middle:end,:), ...
                query(~left),'pchip');
        else
            colors = interp1(knots,base,query,'pchip');
        end
        colors = max(0,min(1,colors));
        method = 'pchip';
    end
    info.provider = 'colorbrewer';
    info.source_palette = scheme.name;
    info.source_url = source.source_url;
    info.authors = source.authors;
    info.license = source.license;
    info.function_name = '';
    info.data_key = scheme.name;
    info.data_file = 'colorbrewer_schemes.json';
    info.palette_type = scheme.type;
    info.sampling_method = method;
    info.matlab_release = '';
end
