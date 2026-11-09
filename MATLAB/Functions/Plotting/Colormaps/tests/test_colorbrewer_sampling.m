function tests = test_colorbrewer_sampling
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    folder = fileparts(fileparts(mfilename('fullpath')));
    old_path = path;
    addpath(folder);
    testCase.addTeardown(@path,old_path);
    testCase.TestData.Source = jsondecode(fileread(fullfile(folder, ...
        'data','colorbrewer_schemes.json')));
end

function testAllDesignedSelections(testCase)
    schemes = testCase.TestData.Source.schemes;
    verifyNumElements(testCase,schemes,35);
    total = 0;
    for ss = 1:numel(schemes)
        for vv = 1:numel(schemes(ss).variants)
            selection = schemes(ss).variants(vv);
            colors = ot_colormap(schemes(ss).name,selection.count, ...
                'Sampling','colorbrewer');
            verifyEqual(testCase,colors,double(selection.colors)/255);
            total = total + 1;
        end
    end
    verifyEqual(testCase,total,265);
end

function testThreeDistinctCategoryColors(testCase)
    expected = [228 26 28;55 126 184;77 175 74]/255;
    [colors,info] = ot_colormap('set1',3,'sampling',"COLORBREWER");
    verifyEqual(testCase,colors,expected);
    verifyEqual(testCase,info.provider,'colorbrewer');
    verifyEqual(testCase,info.source_palette,'Set1');
    verifyEqual(testCase,info.license,'Apache-2.0');
    verifyEqual(testCase,info.sampling_method,'designed');
    verifyNotEqual(testCase,colors,ot_colormap('Set1',3));
end

function testQualitativeRepeatsWithoutBlending(testCase)
    base = ot_colormap('Set1',9,'Sampling','colorbrewer');
    colors = ot_colormap('Set1',12,'Sampling','colorbrewer');
    verifyEqual(testCase,colors,[base;base(1:3,:)]);
    verifyEqual(testCase,ot_colormap('Set1',1,'Sampling','colorbrewer'),base(1,:));
    verifyEqual(testCase,ot_colormap('Set1',2,'Sampling','colorbrewer'),base(1:2,:));
end

function testGradientExtensionsAndSmallCounts(testCase)
    for name = {'Purples','RdBu'}
        largest = 9;
        if strcmp(name{1},'RdBu'), largest = 11; end
        base = ot_colormap(name{1},largest,'Sampling','colorbrewer');
        for N = [32 65 256]
            colors = ot_colormap(name{1},N,'Sampling','colorbrewer');
            verifySize(testCase,colors,[N 3]);
            verifyTrue(testCase,all(isfinite(colors) & colors >= 0 & colors <= 1,'all'));
            verifyEqual(testCase,colors([1 end],:),base([1 end],:));
            if strcmp(name{1},'RdBu') && mod(N,2)
                verifyEqual(testCase,colors((N+1)/2,:),base(6,:));
            end
        end
        small = ot_colormap(name{1},3,'Sampling','colorbrewer');
        verifyEqual(testCase,ot_colormap(name{1},1,'Sampling','colorbrewer'),small(2,:));
        verifyEqual(testCase,ot_colormap(name{1},2,'Sampling','colorbrewer'),small([1 end],:));
        verifySize(testCase,ot_colormap(name{1},0,'Sampling','colorbrewer'),[0 3]);
    end
end

function testDefaultsAndNumericIDs(testCase)
    names = {'Purples','Set1','guppy','parula'};
    for ii = 1:numel(names)
        verifyEqual(testCase,ot_colormap(names{ii},17), ...
            ot_colormap(names{ii},17,'Sampling','default'));
    end
    [~,info] = ot_colormap('Purples',5);
    verifyEqual(testCase,ot_colormap(info.id,5,'Sampling','colorbrewer'), ...
        ot_colormap('Purples',5,'Sampling','colorbrewer'));
end

function testInvalidAndUnsupportedOptions(testCase)
    verifyError(testCase,@() ot_colormap('guppy',5,'Sampling','colorbrewer'), ...
        'ot_colormap:NotColorBrewer');
    verifyError(testCase,@() ot_colormap('parula',5,'Sampling','colorbrewer'), ...
        'ot_colormap:NotColorBrewer');
    verifyError(testCase,@() ot_colormap('Purples',5,'Sampling','other'), ...
        'ot_colormap:InvalidSampling');
    verifyError(testCase,@() ot_colormap('Purples',5,'Other','colorbrewer'), ...
        'ot_colormap:InvalidOptions');
    verifyError(testCase,@() ot_colormap('Purples',5,'Sampling'), ...
        'ot_colormap:InvalidOptions');
end
