function tests = test_ot_colormap
%TEST_OT_COLORMAP Check palette routing, reproducibility and released data.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    plotting_folder = fileparts(fileparts(mfilename('fullpath')));
    previous_path = path;
    addpath(plotting_folder);
    testCase.addTeardown(@path,previous_path);
    testCase.TestData.DataFolder = fullfile(plotting_folder,'data');
    testCase.TestData.Catalog = jsondecode(fileread(fullfile( ...
        testCase.TestData.DataFolder,'palette_catalog.json')));
end

function testNativeMapsUseInstalledMATLAB(testCase)
    % This includes lengths whose native behavior is not equivalent to
    % resampling a stored 256-row table (notably N=0 and N=1).
    names = {'parula','jet','hot','hsv','bone'};
    for nn = 1:numel(names)
        for N = [0 1 64 256]
            expected = feval(names{nn},N);
            actual = ot_colormap(upper(names{nn}),N);
            verifyEqual(testCase,actual,expected);
            verifySize(testCase,actual,[N 3]);
        end
    end
end

function testExternalPaletteAnchorsAndSizes(testCase)
    % Independently recorded endpoints from the audited source tables.
    names = {'guppy','Purples'};
    first_colors = [0.981469 0.567651 0.311148; ...
                    0.988235 0.984314 0.992157];
    last_colors = [0.117753 0.799278 0.750664; ...
                   0.247059 0 0.490196];
    for nn = 1:numel(names)
        for N = [0 1 2 64 256]
            colors = ot_colormap(names{nn},N);
            verifySize(testCase,colors,[N 3]);
            verifyTrue(testCase,isreal(colors));
            verifyTrue(testCase,all(isfinite(colors),'all'));
            verifyTrue(testCase,all(colors >= 0 & colors <= 1,'all'));
            if N >= 2
                verifyEqual(testCase,colors(1,:),first_colors(nn,:),'AbsTol',1e-12);
                verifyEqual(testCase,colors(end,:),last_colors(nn,:),'AbsTol',1e-12);
            elseif N == 1
                % Preserve slanCM's final-endpoint convention for one color.
                verifyEqual(testCase,colors,last_colors(nn,:),'AbsTol',1e-12);
            end
        end
        verifyEqual(testCase,ot_colormap(upper(names{nn}),64), ...
            ot_colormap(names{nn},64));
    end
end

function testDefaultSizeAndOriginalIDs(testCase)
    verifyEqual(testCase,ot_colormap('guppy'),ot_colormap('guppy',256));
    verifySize(testCase,ot_colormap('parula'),[256 3]);
    verifyEqual(testCase,ot_colormap(1,64),ot_colormap('viridis',64));
    verifyEqual(testCase,ot_colormap(6,64),ot_colormap('parula',64));
    verifyEqual(testCase,ot_colormap('  GuPpY  ',64),ot_colormap('guppy',64));
end

function testMetadataIdentifiesActualProvider(testCase)
    [~,native] = ot_colormap('PARULA',16);
    verifyEqual(testCase,native.provider,'matlab');
    verifyEqual(testCase,native.id,6);
    verifyEqual(testCase,native.matlab_release,version('-release'));
    verifyNotEmpty(testCase,native.source_url);
    verifyNotEmpty(testCase,native.license);

    [~,external] = ot_colormap('guppy',16);
    verifyEqual(testCase,external.provider,'external');
    verifyEqual(testCase,external.name,'guppy');
    verifyEqual(testCase,external.id,116);
    verifyEqual(testCase,external.license,'BSD-3-Clause');
    verifyEqual(testCase,external.source_url,'https://github.com/1313e/CMasher');
    verifyEqual(testCase,external.matlab_release,'');
end

function testUnavailableAndUnknownPalettes(testCase)
    verifyError(testCase,@() ot_colormap(192),'ot_colormap:UnknownPaletteID');
    verifyError(testCase,@() ot_colormap('colorblind'),'ot_colormap:UnknownPalette');
    verifyError(testCase,@() ot_colormap(9999),'ot_colormap:UnknownPaletteID');
    verifyError(testCase,@() ot_colormap('not-a-palette'),'ot_colormap:UnknownPalette');
    verifyError(testCase,@() ot_colormap('delete'),'ot_colormap:UnknownPalette');
end

function testInvalidArguments(testCase)
    invalid_sizes = {-1,0.5,NaN,Inf,1+1i,[1 2],[],'64',true};
    for ii = 1:numel(invalid_sizes)
        verifyError(testCase,@() ot_colormap('guppy',invalid_sizes{ii}), ...
            'ot_colormap:InvalidSize');
    end
    invalid_ids = {0,-1,1.5,NaN,Inf,1+1i,[1 2],[]};
    for ii = 1:numel(invalid_ids)
        verifyError(testCase,@() ot_colormap(invalid_ids{ii}), ...
            'ot_colormap:InvalidID');
    end
    invalid_names = {struct(),{'guppy'},["guppy" "parula"],string(missing),true};
    for ii = 1:numel(invalid_names)
        verifyError(testCase,@() ot_colormap(invalid_names{ii}), ...
            'ot_colormap:InvalidName');
    end
end

function testPublishedDataContainsOnlyRegisteredExternalTables(testCase)
    data_file = fullfile(testCase.TestData.DataFolder,'external_palettes.mat');
    variables = whos('-file',data_file);
    verifyEqual(testCase,{variables.name},{'palettes'});
    loaded = load(data_file,'palettes');
    entries = testCase.TestData.Catalog.entries;
    external = strcmp({entries.provider},'external');
    expected_keys = sort(string({entries(external).data_key}).');
    verifyEqual(testCase,sort(string(fieldnames(loaded.palettes))),expected_keys);
    verifyFalse(testCase,any([entries.id] == 192));
    verifyFalse(testCase,any(strcmpi({entries.name},'colorblind')));
    verifyFalse(testCase,isfield(loaded.palettes,'p192'));

    native_entries = entries(strcmp({entries.provider},'matlab'));
    for ii = 1:numel(native_entries)
        key = sprintf('p%03d',native_entries(ii).id);
        verifyFalse(testCase,isfield(loaded.palettes,key));
        verifyEmpty(testCase,native_entries(ii).data_key);
    end
    for ii = find(external)
        colors = loaded.palettes.(entries(ii).data_key);
        verifySize(testCase,colors,[256 3]);
        verifyTrue(testCase,isreal(colors));
        verifyTrue(testCase,all(isfinite(colors),'all'));
        verifyTrue(testCase,all(colors >= 0 & colors <= 1,'all'));
    end
end

function testNativeShadowIsRejected(testCase)
    % Both fixtures restore the path and temporary files after this test.
    folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
    file = fopen(fullfile(folder.Folder,'parula.m'),'w');
    assert(file >= 0,'Could not create the shadow-function fixture.');
    fclose_file = onCleanup(@() fclose(file));
    fprintf(file,'function colors = parula(N)\ncolors = zeros(N,3);\nend\n');
    clear fclose_file
    testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder.Folder));
    verifyError(testCase,@() ot_colormap('parula',16),'ot_colormap:NativeShadowed');
end
