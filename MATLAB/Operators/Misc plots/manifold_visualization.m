clear
clc

%% System definition
folder_path = fullfile(codes_root(),'MATLAB','Data', ...
    'System configurations','Lattice configurations');
file_name = '5_line';

H = lattice_config(false,fullfile(folder_path,file_name));

M = length(H);
sysk = eigk(false,H);
epsilon = sysk.epsilon;

%% Parameters
% Power
P = 3;

% Excitation
exc = 1:3;
N = length(exc);

% Initial conditions
J0 = zeros(M,1);
J0(exc) = P/N;

% Energy
U = epsilon*J0;

% Thermodynamic equilibrium
[~,~,J_eq] = calculate_MC(U,P,epsilon);

%% Plot the manifold!
% J >= 0, sum(J) = P, and epsilon*J = U. Two occupations are
% eliminated outside the displayed axes. Filled circles select up to three
% plotted occupations; sliders fix the additional free occupations.
% Change coordinates edits the selection; Ready applies it.
manifold_figure = plot_manifold(epsilon,P,U,J0,J_eq);

%% Local functions
function fig = plot_manifold(epsilon,P,U,J0,J_eq)
%PLOT_MANIFOLD Explore the nonnegative occupation polytope in three dimensions.
%   Select up to three plotted occupations with the mode circles. Two
%   occupations outside the axes are eliminated using power and energy;
%   sliders fix the remaining coordinates. The complementary pair with
%   the largest energy separation is chosen for each selection. If it has
%   equal energies, exact projection still works but its dimension can drop.
%
%   J0 and J_eq are marked as excitation and equilibrium. J_eq sets the
%   initial slider values. Reset to equilibrium restores J_eq; Reset to
%   excitation restores J0. Degenerate sections are drawn as polygons,
%   segments, or points. If all energies are identical, only one occupation
%   is eliminated because the energy constraint is redundant.
%
%   No Optimization Toolbox is required. FIG.UserData contains mode indices,
%   full and current-section vertices (rows in original modal coordinates),
%   slider bounds, and UI handles for further inspection.

    narginchk(5,5)
    vertices = manifold_vertices(epsilon,P,U);
    if isempty(vertices)
        error('plot_manifold:Infeasible', ...
            'No nonnegative occupations satisfy the specified P and U.')
    end
    epsilon = double(epsilon(:).');
    P = double(P);
    U = double(U);
    M = numel(epsilon);
    eMin = min(epsilon);
    eMax = max(epsilon);
    if eMax == eMin
        constraintCount = 1;
        energies = zeros(1,M);
        energy = 0;
    else
        constraintCount = 2;
        energies = (epsilon-eMin)/(eMax-eMin);
        if P == 0
            energy = 0;
        else
            energy = ((U/P-eMin)/(eMax-eMin))*P;
        end
    end
    axisModes = 1:min(3,M);
    pendingModes = false(1,M);
    pendingModes(axisModes) = true;
    editingCoordinates = false;
    selectionChanged = false;
    savedAxisColors = {};
    savedGraphicColors = cell(0,3);
    savedToolbarVisibility = 'on';
    dependentModes = [];
    sliderModes = [];
    activeModes = [];
    J0 = validate_point(J0,'J0');
    J_eq = validate_point(J_eq,'J_eq');
    fixedJ = J_eq;
    powerScale = P;
    if P == 0
        powerScale = 1;
    end
    nSliders = 0;
    bounds = zeros(0,2);
    sliceVertices = vertices;
    sliceHandle = gobjects(0);
    contextHandle = gobjects(0);
    palette = struct('background',[0.965 0.959 0.943], ...
        'panel',[0.995 0.990 0.974],'plot',[0.987 0.983 0.970], ...
        'ink',[0.16 0.23 0.27],'teal',[0.07 0.59 0.55], ...
        'coral',[0.96 0.38 0.27],'violet',[0.48 0.32 0.74], ...
        'context',[0.64 0.69 0.70]);

    fig = uifigure('Name','Occupation manifold', ...
        'Position',[100 100 1120 720],'Color',palette.background, ...
        'Visible','off','Tag','ManifoldVisualization');
    layout = uigridlayout(fig,[1 2]);
    layout.ColumnWidth = {'1x',320};
    layout.Padding = [12 12 12 12];
    layout.BackgroundColor = palette.background;
    ax = uiaxes(layout);
    ax.Layout.Column = 1;
    hold(ax,'on')
    grid(ax,'on')
    box(ax,'on')
    view(ax,35,25)
    pbaspect(ax,[1 1 1])
    ax.FontSize = 13;
    ax.Color = palette.plot;
    ax.XColor = palette.ink;
    ax.YColor = palette.ink;
    ax.ZColor = palette.ink;
    ax.GridColor = [0.48 0.57 0.59];
    ax.GridAlpha = 0.13;
    ax.Projection = 'perspective';
    % Cool and warm light give the translucent facets a pearly sheen.
    sceneLights = gobjects(2,1);
    sceneLights(1) = light(ax,'Position',[1 -0.4 1.5],'Style','infinite', ...
        'Color',[0.84 0.92 1.00]);
    sceneLights(2) = light(ax,'Position',[-1 0.8 0.7],'Style','infinite', ...
        'Color',[0.60 0.45 0.52]);

    panel = uipanel(layout,'Title','Manifold coordinates','FontSize',16, ...
        'BackgroundColor',palette.panel,'ForegroundColor',palette.ink, ...
        'HighlightColor',[0.84 0.86 0.82]);
    panel.Layout.Column = 2;
    modeColumns = min(5,M);
    modeRows = ceil(M/modeColumns);
    controls = uigridlayout(panel,[9 1]);
    controls.RowHeight = {55,50,min(65*modeRows+8,138),34,42,'1x',22,32,42};
    controls.RowSpacing = 6;
    controls.Padding = [14 10 14 10];
    controls.BackgroundColor = palette.panel;
    uilabel(controls,'Text',sprintf( ...
        '%d modes | %d free coordinates\nP = %.6g     U = %.6g\nJ >= 0', ...
        M,M-constraintCount,P,U),'FontSize',14,'FontColor',palette.ink);
    descriptionLabel = uilabel(controls,'Text','','WordWrap','on','FontColor',palette.ink);
    modePanel = uipanel(controls,'BorderType','none','BackgroundColor',palette.panel);
    modeGrid = uigridlayout(modePanel,[1 1],'Scrollable','on');
    modeGrid.Padding = [0 0 0 0];
    modeGrid.RowHeight = {65*modeRows};
    modeGrid.BackgroundColor = palette.panel;
    modeAxes = uiaxes(modeGrid,'Visible','off','Color','none');
    modeAxes.PositionConstraint = 'outerposition';
    modeAxes.Toolbar.Visible = 'off';
    disableDefaultInteractivity(modeAxes)
    hold(modeAxes,'on')
    xlim(modeAxes,[0.4 modeColumns+0.6])
    ylim(modeAxes,[0.1 modeRows+0.5])
    circles = gobjects(M,1);
    for modeIndex = 1:M
        x = mod(modeIndex-1,modeColumns)+1;
        y = modeRows-floor((modeIndex-1)/modeColumns);
        circles(modeIndex) = scatter(modeAxes,x,y,260, ...
            'MarkerEdgeColor',palette.teal,'LineWidth',1.5, ...
            'MarkerFaceColor','none','HitTest','off','PickableParts','all', ...
            'Tag',sprintf('ManifoldMode%d',modeIndex), ...
            'ButtonDownFcn',@(~,~) toggle_mode(modeIndex));
        text(modeAxes,x,y-0.45,sprintf('J_{%d}',modeIndex), ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',11,'Color',palette.ink,'HitTest','off','PickableParts','none');
    end
    coordinateButton = uibutton(controls,'Text','Change coordinates', ...
        'BackgroundColor',palette.teal,'FontColor',[1 1 1],'FontWeight','bold', ...
        'ButtonPushedFcn',@change_coordinates);
    instructionLabel = uilabel(controls,'Text','','WordWrap','on','FontColor',palette.ink);
    sliderPanel = uipanel(controls,'BorderType','none','BackgroundColor',palette.panel);
    sliderGrid = uigridlayout(sliderPanel,[1 2]);
    sliderGrid.ColumnWidth = {'1x',86};
    sliderGrid.Padding = [6 4 6 4];
    sliderGrid.Scrollable = 'on';
    sliderGrid.BackgroundColor = palette.panel;
    sliders = gobjects(0);
    fields = gobjects(0);
    showProjection = uicheckbox(controls,'Text','Show full projection', ...
        'FontColor',palette.ink, ...
        'Value',M-constraintCount > 3,'ValueChangedFcn',@toggle_projection);
    actions = uigridlayout(controls,[1 2]);
    actions.Padding = [0 0 0 0];
    actions.BackgroundColor = palette.panel;
    resetButton = uibutton(actions,'Text','Reset to equilibrium', ...
        'BackgroundColor',[0.93 0.95 0.91],'FontColor',palette.ink, ...
        'ButtonPushedFcn',@(~,~) set_slice(J_eq));
    excitationButton = uibutton(actions,'Text','Reset to excitation', ...
        'BackgroundColor',[0.93 0.95 0.91],'FontColor',palette.ink, ...
        'ButtonPushedFcn',@(~,~) set_slice(J0));
    status = uilabel(controls,'Text','','WordWrap','on','FontColor',palette.ink);

    initialHandle = plot3(ax,NaN,NaN,NaN,'o', ...
        'MarkerSize',9,'MarkerFaceColor',palette.coral, ...
        'MarkerEdgeColor',[1 0.98 0.94],'LineWidth',1.2);
    initialLabel = text(ax,NaN,NaN,NaN,'excitation','FontSize',11, ...
        'VerticalAlignment','bottom','Clipping','off','Color',palette.ink);
    equilibriumHandle = plot3(ax,NaN,NaN,NaN,'d', ...
        'MarkerSize',9,'MarkerFaceColor',palette.violet, ...
        'MarkerEdgeColor',[1 0.98 0.94],'LineWidth',1.2);
    equilibriumLabel = text(ax,NaN,NaN,NaN,'equilibrium','FontSize',11, ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Clipping','off','Color',palette.ink);
    configure_coordinates()
    fig.Visible = 'on';

    function change_coordinates(~,~)
        if ~editingCoordinates
            pendingModes(:) = false;
            pendingModes(axisModes) = true;
            editingCoordinates = true;
            selectionChanged = false;
            set_plot_editing(true)
            for j = 1:nSliders
                sliders(j).Enable = 'off';
                fields(j).Enable = 'off';
            end
            showProjection.Enable = 'off';
            resetButton.Enable = 'off';
            excitationButton.Enable = 'off';
            instructionLabel.Text = 'Press circles to choose the plotted modes. Filled circles are selected.';
            update_selector()
        elseif selectionChanged && nnz(pendingModes) <= 3
            % Keep the current feasible reference, initially the equilibrium,
            % when a plotted occupation becomes a fixed slider occupation.
            axisModes = find(pendingModes);
            editingCoordinates = false;
            set_plot_editing(false)
            configure_coordinates()
        end
    end

    function toggle_mode(mode)
        if ~editingCoordinates
            return
        end
        pendingModes(mode) = ~pendingModes(mode);
        selectionChanged = true;
        update_selector()
    end

    function update_selector()
        for j = 1:M
            if pendingModes(j)
                circles(j).MarkerFaceColor = palette.teal;
            else
                % A background-colored face keeps the hollow center clickable.
                circles(j).MarkerFaceColor = palette.panel;
            end
            circles(j).HitTest = on_off(editingCoordinates);
        end
        if ~editingCoordinates
            coordinateButton.Text = 'Change coordinates';
            coordinateButton.Enable = 'on';
        elseif selectionChanged && nnz(pendingModes) <= 3
            coordinateButton.Text = 'Ready';
            coordinateButton.Enable = 'on';
        else
            coordinateButton.Text = 'Up to 3 modes may be selected';
            coordinateButton.Enable = 'off';
        end
        data = fig.UserData;
        data.EditingCoordinates = editingCoordinates;
        data.PendingAxisModes = find(pendingModes);
        fig.UserData = data;
    end

    function set_plot_editing(editing)
        colorProperties = {'Color','FaceColor','EdgeColor', ...
            'MarkerFaceColor','MarkerEdgeColor','FaceVertexCData'};
        if editing
            savedAxisColors = get(ax,{'Color','XColor','YColor','ZColor'});
            set(ax,'Color',[0.95 0.95 0.95],'XColor',[0.65 0.65 0.65], ...
                'YColor',[0.65 0.65 0.65],'ZColor',[0.65 0.65 0.65]);
            objects = [sliceHandle;contextHandle;initialHandle;initialLabel; ...
                equilibriumHandle;equilibriumLabel;ax.Title];
            savedGraphicColors = cell(numel(objects)*numel(colorProperties),3);
            colorCount = 0;
            for objectIndex = 1:numel(objects)
                object = objects(objectIndex);
                for propertyIndex = 1:numel(colorProperties)
                    property = colorProperties{propertyIndex};
                    if isprop(object,property)
                        color = get(object,property);
                        if isnumeric(color) && ~isempty(color) && size(color,2) == 3
                            colorCount = colorCount+1;
                            savedGraphicColors(colorCount,:) = {object,property,color};
                            set(object,property,repmat([0.70 0.70 0.70],size(color,1),1));
                        end
                    end
                end
            end
            savedGraphicColors = savedGraphicColors(1:colorCount,:);
            savedToolbarVisibility = ax.Toolbar.Visible;
            ax.Toolbar.Visible = 'off';
            set(sceneLights,'Visible','off')
            disableDefaultInteractivity(ax)
        else
            set(ax,{'Color','XColor','YColor','ZColor'},savedAxisColors);
            for j = 1:size(savedGraphicColors,1)
                set(savedGraphicColors{j,1},savedGraphicColors{j,2},savedGraphicColors{j,3});
            end
            ax.Toolbar.Visible = savedToolbarVisibility;
            set(sceneLights,'Visible','on')
            enableDefaultInteractivity(ax)
        end
    end

    function configure_coordinates()
        remaining = setdiff(1:M,axisModes,'stable');
        if numel(remaining) <= constraintCount
            dependentModes = remaining;
        elseif constraintCount == 1
            dependentModes = remaining(1);
        else
            [~,order] = sort(energies(remaining));
            dependentModes = remaining(order([1 end]));
        end
        sliderModes = setdiff(remaining,dependentModes,'stable');
        activeModes = [axisModes dependentModes];
        nSliders = numel(sliderModes);
        bounds = zeros(nSliders,2);

        labels = {'unused coordinate','unused coordinate','unused coordinate'};
        for k = 1:numel(axisModes)
            labels{k} = sprintf('J_{%d}',axisModes(k));
        end
        xlabel(ax,labels{1}); ylabel(ax,labels{2}); zlabel(ax,labels{3});
        xyz = project(vertices);
        extent = max(max(xyz,[],1),0.08*powerScale);
        xlim(ax,[-0.03*extent(1),1.08*extent(1)])
        ylim(ax,[-0.03*extent(2),1.08*extent(2)])
        zlim(ax,[-0.03*extent(3),1.08*extent(3)])
        pbaspect(ax,[1 1 1])
        point = project(J0);
        initialHandle.XData = point(1);
        initialHandle.YData = point(2);
        initialHandle.ZData = point(3);
        initialLabel.Position = point+0.015*extent;
        equilibriumPoint = project(J_eq);
        equilibriumHandle.XData = equilibriumPoint(1);
        equilibriumHandle.YData = equilibriumPoint(2);
        equilibriumHandle.ZData = equilibriumPoint(3);
        equilibriumLabel.Position = equilibriumPoint-0.015*extent;

        description = sprintf('Axes: %s\nEliminated: %s', ...
            mode_list(axisModes),mode_list(dependentModes));
        if constraintCount == 1
            description = sprintf('%s\nEnergy constraint is redundant.',description);
        elseif numel(dependentModes) == 2 && ...
                abs(diff(energies(dependentModes))) <= 64*eps(1)
            description = sprintf(['Axes: %s\nHidden: %s\n' ...
                'Equal hidden energies; axes are dependent.'], ...
                mode_list(axisModes),mode_list(dependentModes));
        end
        descriptionLabel.Text = description;
        if nSliders == 0
            instructionLabel.Text = 'Full projection onto the displayed axes. No sliders are needed.';
        else
            instructionLabel.Text = 'Sliders fix extra occupations. Each range keeps the section feasible.';
        end
        rebuild_sliders()
        showProjection.Enable = on_off(nSliders > 0);
        resetButton.Enable = on_off(nSliders > 0);
        excitationButton.Enable = on_off(nSliders > 0);
        delete(contextHandle)
        contextHandle = gobjects(0);
        if nSliders > 0
            contextHandle = draw_polytope(ax,xyz,palette.context,0,powerScale);
            contextHandle.LineStyle = '--';
            contextHandle.LineWidth = 0.7;
            if isprop(contextHandle,'FaceColor')
                contextHandle.FaceColor = 'none';
                contextHandle.EdgeColor = palette.context;
            end
            contextHandle.Visible = on_off(showProjection.Value);
        end
        refresh(0)
        update_selector()
    end

    function rebuild_sliders()
        delete(sliderGrid.Children)
        sliders = gobjects(nSliders,1);
        fields = gobjects(nSliders,1);
        if nSliders > 0
            sliderGrid.RowHeight = repmat({27,44},1,nSliders);
        else
            sliderGrid.RowHeight = {'1x'};
        end
        for k = 1:nSliders
            mode = sliderModes(k);
            label = uilabel(sliderGrid,'Text',sprintf('J%d   (energy %.4g)',mode,epsilon(mode)), ...
                'FontColor',palette.ink);
            label.Layout.Row = 2*k-1;
            label.Layout.Column = [1 2];
            sliders(k) = uislider(sliderGrid,'Limits',[0 powerScale], ...
                'Value',fixedJ(mode),'MinorTicks',[],'MajorTicks',[], ...
                'Tag',sprintf('ManifoldSlider%d',mode), ...
                'ValueChangingFcn',@(~,event) change_coordinate(k,event.Value,true), ...
                'ValueChangedFcn',@(~,event) change_coordinate(k,event.Value,false));
            sliders(k).Layout.Row = 2*k;
            sliders(k).Layout.Column = 1;
            fields(k) = uieditfield(sliderGrid,'numeric','Limits',[0 Inf], ...
                'Value',fixedJ(mode),'ValueDisplayFormat','%.5g', ...
                'ValueChangedFcn',@(source,~) change_coordinate(k,source.Value,false));
            fields(k).Layout.Row = 2*k;
            fields(k).Layout.Column = 2;
        end
    end

    function xyz = project(J)
        xyz = zeros(size(J,1),3);
        xyz(:,1:numel(axisModes)) = J(:,axisModes);
    end

    function [J,local] = section(movingMode)
        % Include a moving slider among the unknowns to find its exact range.
        unknown = [activeModes movingMode];
        fixed = setdiff(sliderModes,movingMode,'stable');
        residualP = P-sum(fixedJ(fixed));
        residualU = energy-energies(fixed)*fixedJ(fixed).';
        tolerance = 128*eps(1)*powerScale;
        if residualP < 0 && residualP >= -tolerance
            residualP = 0;
        end
        % Subtracting nearly all the power can leave roundoff larger than
        % the remaining section. Snap only roundoff outside its energy range.
        energyLimits = [min(energies(unknown)),max(energies(unknown))]*residualP;
        if residualU < energyLimits(1) && energyLimits(1)-residualU <= tolerance
            residualU = energyLimits(1);
        elseif residualU > energyLimits(2) && residualU-energyLimits(2) <= tolerance
            residualU = energyLimits(2);
        end
        if residualP < 0
            local = zeros(0,numel(unknown));
        else
            local = manifold_vertices(energies(unknown),residualP,residualU);
        end
        J = repmat(fixedJ,size(local,1),1);
        J(:,unknown) = local;
    end

    function change_coordinate(k,value,dragging)
        if editingCoordinates
            return
        end
        fixedJ(sliderModes(k)) = min(max(value,bounds(k,1)),bounds(k,2));
        refresh(k*dragging)
        % Repair the free occupations after a slider move, retaining the
        % exact slider values and a feasible reference for later selections.
        fixedJ(activeModes) = mean(sliceVertices(:,activeModes),1);
    end

    function set_slice(point)
        if editingCoordinates
            return
        end
        fixedJ = point;
        refresh(0)
    end

    function toggle_projection(~,~)
        if editingCoordinates
            return
        end
        if ~isempty(contextHandle)
            contextHandle.Visible = on_off(showProjection.Value);
        end
    end

    function refresh(draggingIndex)
        [sliceVertices,~] = section([]);
        if isempty(sliceVertices)
            error('plot_manifold:EmptySection','The selected section is infeasible.')
        end
        for j = 1:nSliders
            [~,local] = section(sliderModes(j));
            bounds(j,:) = [min(local(:,end)),max(local(:,end))];
            movable = diff(bounds(j,:)) > 1e-12*powerScale;
            sliders(j).Enable = on_off(movable);
            fields(j).Enable = on_off(movable);
            if movable
                sliders(j).Limits = bounds(j,:);
                sliders(j).MajorTicks = bounds(j,:);
            else
                sliders(j).Limits = [0 powerScale];
                sliders(j).MajorTicks = [];
            end
            if j ~= draggingIndex
                sliders(j).Value = min(max(fixedJ(sliderModes(j)), ...
                    sliders(j).Limits(1)),sliders(j).Limits(2));
            end
            fields(j).Value = fixedJ(sliderModes(j));
        end
        delete(sliceHandle)
        [sliceHandle,dimension] = draw_polytope(ax,project(sliceVertices), ...
            palette.teal,0.45,powerScale);
        if nSliders == 0
            title(ax,sprintf('Occupation polytope (%dD)',dimension))
        else
            title(ax,sprintf('Occupation polytope: %dD section',dimension))
        end
        status.Text = sprintf('Section: %dD, %d vertices\nBoth constraints are satisfied.', ...
            dimension,size(sliceVertices,1));
        fig.UserData = struct('Epsilon',epsilon,'P',P,'U',U,'J0',J0,'J_eq',J_eq, ...
            'DependentModes',dependentModes,'AxisModes',axisModes, ...
            'SliderModes',sliderModes,'Vertices',vertices, ...
            'SliceVertices',sliceVertices,'SliceDimension',dimension, ...
            'SliderValues',fixedJ(sliderModes),'SliderBounds',bounds, ...
            'Axes',ax,'Sliders',sliders,'NumericFields',fields, ...
            'ResetButton',resetButton,'ExcitationButton',excitationButton, ...
            'ProjectionCheckbox',showProjection,'ModeCircles',circles, ...
            'CoordinateButton',coordinateButton,'EditingCoordinates',editingCoordinates, ...
            'PendingAxisModes',find(pendingModes), ...
            'InitialPoint',initialHandle,'InitialLabel',initialLabel, ...
            'EquilibriumPoint',equilibriumHandle,'EquilibriumLabel',equilibriumLabel);
        drawnow limitrate nocallbacks
    end

    function point = validate_point(point,name)
        validateattributes(point,{'numeric'}, ...
            {'vector','real','finite','nonnegative','numel',M},mfilename,name)
        point = double(point(:).');
        if abs(sum(point)-P) > 1e-10*max(P,realmin) || ...
                abs(energies*point.'-energy) > 1e-10*max(P,realmin)
            error('plot_manifold:InvalidOccupationPoint', ...
                '%s must satisfy both the power and energy constraints.',name)
        end
    end
end

function [handle,dimension] = draw_polytope(ax,xyz,color,alpha,scale)
    xyz = uniquetol(xyz/scale,1e-12,'ByRows',true,'DataScale',1)*scale;
    centered = (xyz-mean(xyz,1))/scale;
    [~,S,basis] = svd(centered,0);
    if size(xyz,1) == 1
        dimension = 0;
    else
        dimension = sum(diag(S) > 1e-11);
    end
    switch dimension
        case 0
            handle = plot3(ax,xyz(1,1),xyz(1,2),xyz(1,3),'.', ...
                'Color',color,'MarkerSize',22);
        case 1
            coordinate = centered*basis(:,1);
            [~,lo] = min(coordinate);
            [~,hi] = max(coordinate);
            ends = xyz([lo hi],:);
            handle = plot3(ax,ends(:,1),ends(:,2),ends(:,3),'-o', ...
                'Color',color,'LineWidth',2,'MarkerFaceColor',color);
        case 2
            plane = centered*basis(:,1:2);
            faces = convhull(plane(:,1),plane(:,2)).';
            handle = patch(ax,'Vertices',xyz,'Faces',faces, ...
                'FaceColor',[0.80 0.86 0.96],'FaceAlpha',alpha, ...
                'EdgeColor',[0.24 0.40 0.51],'EdgeAlpha',0.75,'LineWidth',1.2);
        otherwise
            [faces,normals] = polytope_faces(xyz/scale);
            % Soft orientation-dependent tints suggest opal while preserving
            % flat faces. Transparency keeps interior points and edges visible.
            tints = [0.80 0.83 1.00; 0.68 0.91 0.91; 1.00 0.84 0.90];
            weights = abs(normals).^1.5;
            weights = weights./sum(weights,2);
            faceColors = 0.12+0.88*(weights*tints);
            handle = patch(ax,'Vertices',xyz,'Faces',faces, ...
                'FaceColor','flat','FaceVertexCData',faceColors,'FaceAlpha',alpha, ...
                'EdgeColor',[0.24 0.40 0.51],'EdgeAlpha',0.75, ...
                'LineWidth',1.2,'FaceLighting','flat','EdgeLighting','none', ...
                'AmbientStrength',0.70,'DiffuseStrength',0.45, ...
                'SpecularStrength',0.75,'SpecularExponent',22);
    end
end

function [faces,facetNormals] = polytope_faces(xyz)
% Merge coplanar hull triangles so drawn edges are actual polytope edges.
    triangles = convhull(xyz(:,1),xyz(:,2),xyz(:,3));
    edge1 = xyz(triangles(:,2),:)-xyz(triangles(:,1),:);
    edge2 = xyz(triangles(:,3),:)-xyz(triangles(:,1),:);
    normals = cross(edge1,edge2,2);
    normals = normals./vecnorm(normals,2,2);
    centers = (xyz(triangles(:,1),:)+xyz(triangles(:,2),:)+xyz(triangles(:,3),:))/3;
    inward = sum(normals.*(centers-mean(xyz,1)),2) < 0;
    normals(inward,:) = -normals(inward,:);
    offsets = sum(normals.*centers,2);
    remaining = true(size(triangles,1),1);
    polygons = cell(size(triangles,1),1);
    facetNormals = zeros(size(triangles,1),3);
    count = 0;
    while any(remaining)
        first = find(remaining,1);
        coplanar = remaining & vecnorm(normals-normals(first,:),2,2) < 1e-9 & ...
            abs(offsets-offsets(first)) < 1e-10;
        indices = unique(triangles(coplanar,:));
        tangent = edge1(first,:)/norm(edge1(first,:));
        plane = xyz(indices,:)*[tangent;cross(normals(first,:),tangent)].';
        boundary = convhull(plane(:,1),plane(:,2));
        count = count+1;
        facetNormals(count,:) = normals(first,:);
        polygons{count} = indices(boundary(1:end-1)).';
        remaining(coplanar) = false;
    end
    polygons = polygons(1:count);
    facetNormals = facetNormals(1:count,:);
    faces = nan(count,max(cellfun(@numel,polygons)));
    for k = 1:count
        faces(k,1:numel(polygons{k})) = polygons{k};
    end
end

function text = mode_list(modes)
    if isempty(modes)
        text = '(none)';
    else
        text = strjoin(arrayfun(@(k) sprintf('J%d',k),modes, ...
            'UniformOutput',false),', ');
    end
end

function value = on_off(condition)
    if condition
        value = 'on';
    else
        value = 'off';
    end
end

function vertices = manifold_vertices(epsilon, P, U)
%MANIFOLD_VERTICES Vertices of a nonnegative modal-occupation polytope.
%   V = MANIFOLD_VERTICES(EPSILON,P,U) returns one vertex per row of
%       {J >= 0 : sum(J) = P, EPSILON(:).' * J = U}.
%   An infeasible system returns a 0-by-numel(EPSILON) matrix. No
%   Optimization Toolbox is required. EPSILON may be a row or column vector.
%
%   With two independent equalities every vertex has at most two nonzero
%   occupations. Enumerating pairs is therefore exact, including repeated
%   eigenvalues. Energies are shifted and scaled before comparisons.

    narginchk(3, 3);
    validateattributes(epsilon, {'numeric'}, ...
        {'real', 'finite', 'vector', 'nonempty'}, mfilename, 'epsilon', 1);
    validateattributes(P, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'P', 2);
    validateattributes(U, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'U', 3);

    epsilon = double(epsilon(:).');
    P = double(P);
    U = double(U);
    modeCount = numel(epsilon);
    vertices = zeros(0, modeCount);

    if P == 0
        if U == 0
            vertices = zeros(1, modeCount);
        end
        return
    end

    meanEnergy = U / P;
    if ~isfinite(meanEnergy)
        return
    end
    energyMin = min(epsilon);
    energySpan = max(epsilon) - energyMin;
    tolerance = 64 * eps(1);

    if energySpan == 0
        % Here the energy constraint is redundant or inconsistent. Avoid an
        % absolute tolerance floor: zero energies cannot support nonzero U.
        energyScale = max(abs([meanEnergy, energyMin]));
        if abs(meanEnergy - energyMin) <= tolerance * energyScale
            vertices = P * eye(modeCount);
        end
        return
    end

    if isinf(energySpan)
        % Finite endpoints of opposite sign can have an overflowing span.
        % Scale first in this exceptional case, then shift normally.
        energyScale = max(abs(epsilon));
        epsilon = epsilon / energyScale;
        meanEnergy = meanEnergy / energyScale;
        energyMin = min(epsilon);
        energySpan = max(epsilon) - energyMin;
    end
    normalizedEnergy = (epsilon - energyMin) / energySpan;
    target = (meanEnergy - energyMin) / energySpan;
    if ~isfinite(target) || target < -tolerance || target > 1 + tolerance
        return
    end
    target = min(1, max(0, target));

    % Near-target modes represent singleton vertices. Do not also generate
    % pairs involving these modes, which would duplicate their endpoints.
    onTarget = find(abs(normalizedEnergy - target) <= tolerance);
    below = find(normalizedEnergy < target - tolerance);
    above = find(normalizedEnergy > target + tolerance);
    vertexCount = numel(onTarget) + numel(below) * numel(above);
    vertices = zeros(vertexCount, modeCount);
    row = 0;
    for mode = onTarget
        row = row + 1;
        vertices(row, mode) = P;
    end
    for lowMode = below
        for highMode = above
            row = row + 1;
            fraction = (normalizedEnergy(highMode) - target) / ...
                (normalizedEnergy(highMode) - normalizedEnergy(lowMode));
            vertices(row, lowMode) = P * fraction;
            vertices(row, highMode) = P - vertices(row, lowMode);
        end
    end
end
