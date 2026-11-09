function [H_lattice,B] = lattice_config(plotting,load_location)
%LATTICE_CONFIG Assemble a lattice from independently editable sublattices.
%   H_LATTICE = LATTICE_CONFIG() opens a blocking graphical editor and
%   returns only the full composite Hamiltonian.
%   Add creates a fresh sublattice, Edit reopens the selected sublattice, and
%   Delete removes it. Sublattices translate as rigid bodies in the outer
%   editing plane and rotate about their local origins; their local
%   coordinates and internal Hamiltonians remain unchanged.
%
%   [H_LATTICE,B] = LATTICE_CONFIG() additionally returns the ordered
%   site indices of each subsystem. B{k} selects Subsystem k from the packed
%   Hamiltonian as H_LATTICE(B{k},B{k}). H_LATTICE is ordered with every site
%   in Subsystem 1 first, followed by every site in Subsystem 2, and so on.
%   Cancel or closing the editor returns H_LATTICE = [] and B = {}.
%
%   H_LATTICE = LATTICE_CONFIG(PLOTTING) additionally opens a
%   maximized interactive modal viewer after a successful Build when the
%   logical scalar PLOTTING is true. PLOTTING defaults to false. The viewer
%   switches among the full lattice and its ordered sublattices, and among
%   the eigenvalue spectrum, modal intensities, and modal phases.
%
%   H_LATTICE = LATTICE_CONFIG(LOAD_LOCATION) loads a saved composite
%   lattice without opening the editor. Requesting a second output also
%   returns its subsystem index sets. LOAD_LOCATION may be a bare file name,
%   relative path, or absolute path. A .mat extension is appended when
%   omitted. LATTICE_CONFIG(PLOTTING,LOAD_LOCATION) also opens the
%   modal viewer after loading when PLOTTING is true.
%
%   INTER-SUBLATTICE COUPLINGS
%   Each sublattice has its own internal shortest edge and coupling law.
%   Neither is changed by this editor. Couplings between different
%   sublattices instead use the outer profile, absolute range cutoff, and
%   decay parameter. A cross-sublattice pair is coupled only when its
%   separation does not exceed that absolute cutoff.
%   The outer coupling law is normalized by
%
%       d_min = min |x_a - x_b|,   a ~= b,
%
%   over every distinct pair of sites in the assembled lattice, including
%   pairs within the same sublattice. Here x_a and x_b are global positions
%   after sublattice translation and rotation. The outer d_min can therefore
%   change when a sublattice moves, but internal Hamiltonian blocks retain
%   their own independently configured coupling laws. With q = d/d_min, the
%   available outer profiles are
%
%       Power law:        kappa = q^(-gamma)
%       Exponential:      kappa = exp[-gamma(q-1)]
%       Fiber asymptotic: kappa = q^(-1/2) exp[-gamma(q-1)].
%
%   Cross-sublattice couplings are real and reciprocal. Internal blocks are
%   inserted into H_LATTICE without modification, so internal phases,
%   asymmetry, manual couplings, detunings, and disorder are all preserved.
%
%   SAVED LATTICES
%   Save records the complete editable composite model in a MAT-file under
%   the variable name latticeConfiguration. Load restores sublattice models,
%   rigid transforms, and outer coupling properties before rebuilding the
%   composite Hamiltonian. Supplying LOAD_LOCATION performs this same rebuild
%   without constructing either editor.
%
%   MODAL VIEWER
%   When PLOTTING is true, each selected Hamiltonian is diagonalized only on
%   first use and then cached. Modes are ordered lexicographically by the real
%   and imaginary parts of their eigenvalues. Spectrum tiles show both parts;
%   intensity tiles color sites by abs(psi).^2 using a hot map truncated at
%   95 percent; phase tiles use angle(psi) and hsv. Profile tiles retain the
%   selected lattice geometry and neutral versions of edited edge curves,
%   dotted segments, and directional path separation.

    %% Inputs, outputs, and composite state

    narginchk(0, 2);

    % Supported positional forms are (), (plotting), (load_location), and
    % (plotting,load_location). A single text input is therefore dispatched
    % as a load request rather than as the plotting flag.
    loadRequested = false;
    if nargin == 0
        plotting = false;
    elseif nargin == 1 && islogical(plotting) && isscalar(plotting)
        % A lone logical controls plotting for an interactive build.
    elseif nargin == 1 && (ischar(plotting) ...
            || (isstring(plotting) && isscalar(plotting)))
        % Treat a single text input as LOAD_LOCATION.
        load_location = plotting;
        plotting = false;
        loadRequested = true;
    elseif nargin == 1
        error('lattice_config:InvalidPlottingOption', ...
            'Plotting must be true or false.');
    else
        if ~islogical(plotting) || ~isscalar(plotting)
            error('lattice_config:InvalidPlottingOption', ...
                'Plotting must be true or false.');
        end
        loadRequested = true;
    end

    H_lattice = [];
    B = cell(0, 1);
    internalBlocks = cell(0, 1);
    configurationAccepted = false;

    regularFontSize = 14;
    panelTitleFontSize = 22;
    buttonFontSize = 16;
    axisLabelFontSize = 14;
    axisTitleFontSize = 22;
    sectionHeaderHeight = 38;
    sectionSeparatorColor = [0.55 0.55 0.55];

    % A composite record owns rigid placement plus a complete embedded model
    % and Hamiltonian for each sublattice. Outer edits never mutate local site
    % coordinates or internal couplings; Edit replaces the embedded pair as a
    % unit after the child editor successfully returns.
    state = struct();
    state.sublattices = struct( ...
        'uid', {}, ...
        'origin', {}, ...
        'rotationDegrees', {}, ...
        'configuration', {}, ...
        'H', {});
    state.nextUid = 1;
    state.selectedUid = [];
    state.gamma = 1;
    state.cutoffDistance = Inf;
    state.couplingProfile = "Exponential";
    state.globalDMin = nan;
    state.interaction = defaultInteractionState();
    state.initialized = false;
    state.cancelRequested = false;
    defaultPersistenceDirectory = string(fullfile(codes_root(), ...
        'MATLAB','Data','System configurations','Lattice configurations'));
    defaultPersistenceFileName = "lattice";
    state.persistenceMode = "path";
    state.pendingPersistenceAction = "";

    % persistenceMode drives a small UI state machine:
    % path -> confirm -> success/failure -> path.

    % Undo snapshots contain semantic model state only. Selection, camera
    % position, and persistence-interface status remain outside the history.
    undoStack = cell(0, 1);
    undoLimit = 40;
    dragUndoSnapshot = [];
    dragUndoRecorded = false;

    if loadRequested
        savedPath = resolveSavedLatticePath(load_location);
        loadedConfiguration = readLatticeConfigurationFile(savedPath);
        installConfigurationModel(loadedConfiguration);
        [H_lattice,B,internalBlocks] = calculateCompositeHamiltonian();
        if plotting
            launchModalViewer();
        end
        return;
    end

    %% Main window and column layout

    leftColumnWidth = 398;
    leftRowHeights = [161 44 162 415 163];
    leftRowSpacing = 5;
    mainPadding = [6 6 6 6];
    mainColumnSpacing = 6;
    defaultFigurePosition = defaultEditorWindowPosition( ...
        [141 1032], leftColumnWidth, leftRowHeights, leftRowSpacing, ...
        mainPadding, mainColumnSpacing);

    ui = struct();
    ui.figure = uifigure( ...
        'Name', 'DNLSE lattice configuration', ...
        'Tag', 'LatticeConfigurationFigure', ...
        'Position', defaultFigurePosition, ...
        'Color', [0.96 0.97 0.98], ...
        'Visible', 'off', ...
        'WindowButtonMotionFcn', @pointerMoved, ...
        'WindowButtonUpFcn', @pointerReleased, ...
        'WindowScrollWheelFcn', @zoomEditingPlane, ...
        'WindowKeyPressFcn', @handleKeyPress, ...
        'CloseRequestFcn', @cancelConfiguration);

    mainGrid = uigridlayout(ui.figure, [1 2]);
    mainGrid.ColumnWidth = {leftColumnWidth, '1x'};
    mainGrid.RowHeight = {'1x'};
    mainGrid.Padding = mainPadding;
    mainGrid.ColumnSpacing = mainColumnSpacing;

    leftGrid = uigridlayout(mainGrid, [5 1]);
    leftGrid.Layout.Row = 1;
    leftGrid.Layout.Column = 1;
    leftGrid.RowHeight = num2cell(leftRowHeights);
    leftGrid.ColumnWidth = {'1x'};
    leftGrid.Padding = [0 0 0 0];
    leftGrid.RowSpacing = leftRowSpacing;
    leftGrid.Scrollable = 'on';

    %% Inter-sublattice coupling properties

    [latticePanel,latticeBody,latticeTitleLabel] = createSectionPanel( ...
        leftGrid, 'Lattice properties', 'LatticePropertiesPanel');
    latticePanel.Layout.Row = 1;

    latticeGrid = uigridlayout(latticeBody, [4 1]);
    latticeGrid.ColumnWidth = {'1x'};
    latticeGrid.RowHeight = {30, 25, 26, 26};
    latticeGrid.Padding = [7 2 7 3];
    latticeGrid.RowSpacing = 3;
    latticeGrid.ColumnSpacing = 0;

    profileRowGrid = uigridlayout(latticeGrid, [1 2]);
    profileRowGrid.Layout.Row = 2;
    profileRowGrid.ColumnWidth = {'1x', '1x'};
    profileRowGrid.RowHeight = {'1x'};
    profileRowGrid.Padding = [0 0 0 0];
    profileRowGrid.ColumnSpacing = 5;
    profileLabel = uilabel(profileRowGrid, 'Text', 'Coupling profile');
    profileLabel.Layout.Column = 1;
    ui.profileDropdown = uidropdown(profileRowGrid, ...
        'Items', {'Power law', 'Exponential', 'Fiber asymptotic'}, ...
        'Value', char(state.couplingProfile), ...
        'Tag', 'CouplingProfileDropdown', ...
        'ValueChangedFcn', @profileChanged);
    ui.profileDropdown.Layout.Column = 2;

    ui.formulaPanel = uipanel(latticeGrid, ...
        'BorderType', 'line', ...
        'BorderWidth', 1, ...
        'BorderColor', [0.72 0.72 0.72], ...
        'BackgroundColor', [0.97 0.975 0.98], ...
        'Tag', 'InterSublatticeLawPanel');
    ui.formulaPanel.Layout.Row = 1;
    formulaGrid = uigridlayout(ui.formulaPanel, [1 1]);
    formulaGrid.RowHeight = {'1x'};
    formulaGrid.ColumnWidth = {'1x'};
    formulaGrid.Padding = [0 0 0 0];
    ui.formulaLabel = uilabel(formulaGrid, ...
        'Text', couplingFormula(), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.97 0.975 0.98], ...
        'FontColor', [0.22 0.24 0.28], ...
        'Tag', 'InterSublatticeLawLabel');

    cutoffRowGrid = uigridlayout(latticeGrid, [1 2]);
    cutoffRowGrid.Layout.Row = 3;
    cutoffRowGrid.ColumnWidth = {'1x', '1x'};
    cutoffRowGrid.RowHeight = {'1x'};
    cutoffRowGrid.Padding = [0 0 0 0];
    cutoffRowGrid.ColumnSpacing = 5;
    ui.cutoffLabel = uilabel(cutoffRowGrid, ...
        'Text', cutoffLabelText(), ...
        'Tooltip', ['This absolute cutoff governs only cross-sublattice ' ...
                    'couplings. Their normalization uses the globally ' ...
                    'shortest pair distance.']);
    ui.cutoffLabel.Layout.Column = 1;
    ui.cutoffField = uieditfield(cutoffRowGrid, 'numeric', ...
        'Value', state.cutoffDistance, ...
        'Limits', [0 Inf], ...
        'Tag', 'InterSublatticeCutoffField', ...
        'Tooltip', ['Maximum absolute separation for couplings between ' ...
                    'different sublattices. Values entered below the ' ...
                    'current global d_min are raised to d_min.'], ...
        'ValueChangedFcn', @cutoffChanged);
    ui.cutoffField.Layout.Column = 2;

    gammaRowGrid = uigridlayout(latticeGrid, [1 2]);
    gammaRowGrid.Layout.Row = 4;
    gammaRowGrid.ColumnWidth = {'1x', '1x'};
    gammaRowGrid.RowHeight = {'1x'};
    gammaRowGrid.Padding = [0 0 0 0];
    gammaRowGrid.ColumnSpacing = 5;
    gammaLabel = uilabel(gammaRowGrid, 'Text', 'Decay parameter, γ');
    gammaLabel.Layout.Column = 1;
    ui.gammaField = uieditfield(gammaRowGrid, 'numeric', ...
        'Value', state.gamma, ...
        'Limits', [eps Inf], ...
        'Tag', 'InterSublatticeGammaField', ...
        'Tooltip', 'Distance decay for couplings between sublattices', ...
        'ValueChangedFcn', @gammaChanged);
    ui.gammaField.Layout.Column = 2;

    %% Sublattice creation and selected-sublattice properties

    ui.addButton = uibutton(leftGrid, 'push', ...
        'Text', 'Add sublattice', ...
        'Tag', 'AddSublatticeButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.84 0.96 0.98], ...
        'Tooltip', 'Create a sublattice at the lattice origin', ...
        'ButtonPushedFcn', @addSublattice);
    ui.addButton.Layout.Row = 2;

    [ui.sublatticePanel,sublatticeBody,ui.sublatticeTitleLabel] = ...
        createSectionPanel(leftGrid, 'Sublattice properties', ...
        'SelectedSublatticePropertiesPanel');
    ui.sublatticePanel.Layout.Row = 3;

    sublatticeGrid = uigridlayout(sublatticeBody, [4 2]);
    sublatticeGrid.ColumnWidth = {'1x', 100};
    sublatticeGrid.RowHeight = {26, 26, 26, 30};
    sublatticeGrid.Padding = [7 2 7 3];
    sublatticeGrid.RowSpacing = 3;
    sublatticeGrid.ColumnSpacing = 5;

    xLabel = uilabel(sublatticeGrid, 'Text', 'X-coordinate');
    xLabel.Layout.Row = 1;
    xLabel.Layout.Column = 1;
    ui.xField = uieditfield(sublatticeGrid, 'numeric', ...
        'Value', 0, ...
        'Enable', 'off', ...
        'Tag', 'SublatticeXField', ...
        'Tooltip', 'X-coordinate of the selected sublattice origin', ...
        'ValueChangedFcn', @(source,event) originFieldChanged(source,event,1));
    ui.xField.Layout.Row = 1;
    ui.xField.Layout.Column = 2;

    yLabel = uilabel(sublatticeGrid, 'Text', 'Y-coordinate');
    yLabel.Layout.Row = 2;
    yLabel.Layout.Column = 1;
    ui.yField = uieditfield(sublatticeGrid, 'numeric', ...
        'Value', 0, ...
        'Enable', 'off', ...
        'Tag', 'SublatticeYField', ...
        'Tooltip', 'Y-coordinate of the selected sublattice origin', ...
        'ValueChangedFcn', @(source,event) originFieldChanged(source,event,2));
    ui.yField.Layout.Row = 2;
    ui.yField.Layout.Column = 2;

    rotationLabel = uilabel(sublatticeGrid, 'Text', 'Rotation (°)');
    rotationLabel.Layout.Row = 3;
    rotationLabel.Layout.Column = 1;
    ui.rotationField = uieditfield(sublatticeGrid, 'numeric', ...
        'Value', 0, ...
        'Enable', 'off', ...
        'Tag', 'SublatticeRotationField', ...
        'Tooltip', 'Counterclockwise rotation about the selected sublattice origin', ...
        'ValueChangedFcn', @rotationFieldChanged);
    ui.rotationField.Layout.Row = 3;
    ui.rotationField.Layout.Column = 2;

    ui.sublatticeHint = uilabel(sublatticeGrid, ...
        'Text', 'Select a sublattice in the editing plane', ...
        'WordWrap', 'on', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'center', ...
        'FontColor', [0.42 0.45 0.50]);
    ui.sublatticeHint.Layout.Row = [1 4];
    ui.sublatticeHint.Layout.Column = [1 2];

    ui.sublatticeActionGrid = uigridlayout(sublatticeGrid, [1 2]);
    ui.sublatticeActionGrid.Layout.Row = 4;
    ui.sublatticeActionGrid.Layout.Column = [1 2];
    ui.sublatticeActionGrid.ColumnWidth = {'1x', '1x'};
    ui.sublatticeActionGrid.RowHeight = {'1x'};
    ui.sublatticeActionGrid.Padding = [0 0 0 0];
    ui.sublatticeActionGrid.ColumnSpacing = 5;
    ui.deleteButton = uibutton(ui.sublatticeActionGrid, 'push', ...
        'Text', 'Delete', ...
        'Enable', 'off', ...
        'Tag', 'DeleteSublatticeButton', ...
        'BackgroundColor', [1.00 0.86 0.86], ...
        'Tooltip', 'Delete the selected sublattice', ...
        'ButtonPushedFcn', @deleteSelectedSublattice);
    ui.deleteButton.Layout.Column = 1;
    ui.editButton = uibutton(ui.sublatticeActionGrid, 'push', ...
        'Text', 'Edit', ...
        'Enable', 'off', ...
        'Tag', 'EditSublatticeButton', ...
        'BackgroundColor', [0.84 0.95 0.84], ...
        'Tooltip', 'Edit the selected sublattice', ...
        'ButtonPushedFcn', @editSublattice);
    ui.editButton.Layout.Column = 2;
    ui.sublatticeSelectionControls = {xLabel, ui.xField, yLabel, ui.yField, ...
        rotationLabel, ui.rotationField, ui.sublatticeActionGrid};

    %% Instructions, build actions, and persistence UI

    [controlsPanel,controlsBody,controlsTitleLabel] = createSectionPanel( ...
        leftGrid, 'Controls', 'LatticeControlsPanel');
    controlsPanel.Layout.Row = 4;
    controlsGrid = uigridlayout(controlsBody, [6 1]);
    controlsGrid.RowHeight = {22, '1x', 22, '1x', 22, '1x'};
    controlsGrid.ColumnWidth = {'1x'};
    controlsGrid.Padding = [8 6 8 6];
    controlsGrid.RowSpacing = 4;
    sublatticeInstantiationHeading = uilabel(controlsGrid, ...
        'Text', 'Sublattice instantiation', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
    sublatticeInstantiationHeading.Layout.Row = 1;
    sublatticeInstantiationDescription = uihtml(controlsGrid, ...
        'HTMLSource', justifiedControlsHtml( ...
            ['‘Add sublattice’ opens an external sublattice configuration ' ...
             'window that can be used to build a sublattice at the origin. ' ...
             '‘Edit’ and ‘Delete’ can then be used to edit or delete it. ' ...
             'Moving, rotating, or changing outer coupling properties does ' ...
             'not modify a sublattice’s internal model.']), ...
        'Tag', 'LatticeControlsLabel');
    sublatticeInstantiationDescription.Layout.Row = 2;
    editingPlaneHeading = uilabel(controlsGrid, ...
        'Text', 'Editing plane', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
    editingPlaneHeading.Layout.Row = 3;
    editingPlaneDescription = uihtml(controlsGrid, ...
        'HTMLSource', justifiedControlsHtml( ...
            ['‘Left-click’ to select a sublattice and drag to move it. ' ...
             '‘Delete’ or ‘Backspace’ can be used to remove sublattices, ' ...
             'and ‘Ctrl+Z’ undoes the previous lattice change. In the ' ...
             'editing plane, ‘Scroll’ to zoom, use the ‘Hand icon’ to pan, ' ...
             'and the ‘House icon’ to recenter.']), ...
        'Tag', 'LatticeEditingControlsLabel');
    editingPlaneDescription.Layout.Row = 4;
    couplingPropertiesHeading = uilabel(controlsGrid, ...
        'Text', 'Coupling properties', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
    couplingPropertiesHeading.Layout.Row = 5;
    couplingPropertiesDescription = uihtml(controlsGrid, ...
        'HTMLSource', justifiedControlsHtml( ...
            ['Modifying the couplings between whole sublattices is ' ...
             'accomplished by adjusting the coupling profile, range cutoff, ' ...
             'and decay parameter. The minimum distance parameter is set by ' ...
             'the closest pair of nodes, globally.']), ...
        'Tag', 'LatticePropertiesControlsLabel');
    couplingPropertiesDescription.Layout.Row = 6;

    actionGrid = uigridlayout(leftGrid, [3 1]);
    actionGrid.Layout.Row = 5;
    actionGrid.RowHeight = {44, 44, 65};
    actionGrid.ColumnWidth = {'1x'};
    actionGrid.Padding = [0 0 0 0];
    actionGrid.RowSpacing = 5;

    buildButtonGrid = uigridlayout(actionGrid, [1 2]);
    buildButtonGrid.Layout.Row = 1;
    buildButtonGrid.ColumnWidth = {'1x', '1x'};
    buildButtonGrid.RowHeight = {'1x'};
    buildButtonGrid.Padding = [0 0 0 0];
    buildButtonGrid.ColumnSpacing = 8;
    ui.cancelButton = uibutton(buildButtonGrid, 'push', ...
        'Text', 'Cancel', ...
        'Tag', 'CancelLatticeButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [1.00 0.86 0.86], ...
        'ButtonPushedFcn', @cancelConfiguration);
    ui.acceptButton = uibutton(buildButtonGrid, 'push', ...
        'Text', 'Build', ...
        'Tag', 'BuildLatticeButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.84 0.95 0.84], ...
        'ButtonPushedFcn', @acceptConfiguration);

    persistenceButtonGrid = uigridlayout(actionGrid, [1 2]);
    persistenceButtonGrid.Layout.Row = 2;
    persistenceButtonGrid.ColumnWidth = {'1x', '1x'};
    persistenceButtonGrid.RowHeight = {'1x'};
    persistenceButtonGrid.Padding = [0 0 0 0];
    persistenceButtonGrid.ColumnSpacing = 8;
    ui.loadButton = uibutton(persistenceButtonGrid, 'push', ...
        'Text', 'Load', ...
        'Tag', 'LoadLatticeButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.98 0.96 1.00], ...
        'ButtonPushedFcn', @(~,~) beginPersistenceAction("load"));
    ui.saveButton = uibutton(persistenceButtonGrid, 'push', ...
        'Text', 'Save', ...
        'Tag', 'SaveLatticeButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.96 0.98 1.00], ...
        'ButtonPushedFcn', @(~,~) beginPersistenceAction("save"));

    confirmationButtonGrid = uigridlayout(actionGrid, [1 2]);
    confirmationButtonGrid.Layout.Row = 2;
    confirmationButtonGrid.ColumnWidth = {'1x', 44};
    confirmationButtonGrid.RowHeight = {'1x'};
    confirmationButtonGrid.Padding = [0 0 0 0];
    confirmationButtonGrid.ColumnSpacing = 8;
    confirmationButtonGrid.Visible = 'off';
    ui.confirmDirectoryButton = uibutton(confirmationButtonGrid, 'push', ...
        'Text', 'Confirm directory?', ...
        'Tag', 'ConfirmLatticeDirectoryButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.55 0.80 0.55], ...
        'ButtonPushedFcn', @confirmPersistenceAction);
    ui.confirmDirectoryButton.Layout.Column = 1;
    ui.cancelDirectoryButton = uibutton(confirmationButtonGrid, 'push', ...
        'Text', 'X', ...
        'Tag', 'CancelLatticeDirectoryButton', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.90 0.48 0.48], ...
        'ButtonPushedFcn', @cancelPersistenceAction);
    ui.cancelDirectoryButton.Layout.Column = 2;

    persistenceAreaGrid = uigridlayout(actionGrid, [3 1]);
    persistenceAreaGrid.Layout.Row = 3;
    persistenceAreaGrid.RowHeight = {'1x', 30, '1x'};
    persistenceAreaGrid.ColumnWidth = {'1x'};
    persistenceAreaGrid.Padding = [0 0 0 0];
    persistenceAreaGrid.RowSpacing = 0;

    persistenceInputGrid = uigridlayout(persistenceAreaGrid, [2 2]);
    persistenceInputGrid.Layout.Row = [1 3];
    persistenceInputGrid.ColumnWidth = {88, '1x'};
    persistenceInputGrid.RowHeight = {30, 30};
    persistenceInputGrid.Padding = [0 0 0 0];
    persistenceInputGrid.RowSpacing = 5;
    persistenceInputGrid.ColumnSpacing = 5;
    folderLabel = uilabel(persistenceInputGrid, 'Text', 'Folder path');
    folderLabel.Layout.Row = 1;
    folderLabel.Layout.Column = 1;
    ui.folderPathField = uieditfield(persistenceInputGrid, 'text', ...
        'Value', char(defaultPersistenceDirectory), ...
        'Tag', 'LatticeFolderPathField', ...
        'HorizontalAlignment', 'left', ...
        'Tooltip', 'Folder containing the composite-lattice MAT-file', ...
        'ValueChangedFcn', @persistenceLocationChanged);
    ui.folderPathField.Layout.Row = 1;
    ui.folderPathField.Layout.Column = 2;
    fileLabel = uilabel(persistenceInputGrid, 'Text', 'File name');
    fileLabel.Layout.Row = 2;
    fileLabel.Layout.Column = 1;
    ui.fileNameField = uieditfield(persistenceInputGrid, 'text', ...
        'Value', char(defaultPersistenceFileName), ...
        'Tag', 'LatticeFileNameField', ...
        'HorizontalAlignment', 'left', ...
        'Tooltip', 'The .mat extension is added automatically when omitted', ...
        'ValueChangedFcn', @persistenceLocationChanged);
    ui.fileNameField.Layout.Row = 2;
    ui.fileNameField.Layout.Column = 2;

    persistenceStatusGrid = uigridlayout(persistenceAreaGrid, [1 2]);
    persistenceStatusGrid.Layout.Row = 2;
    persistenceStatusGrid.ColumnWidth = {'1x', 30};
    persistenceStatusGrid.RowHeight = {30};
    persistenceStatusGrid.Padding = [0 0 0 0];
    persistenceStatusGrid.ColumnSpacing = 5;
    persistenceStatusGrid.Visible = 'off';
    ui.persistenceStatusField = uieditfield(persistenceStatusGrid, 'text', ...
        'Value', 'Saved', ...
        'Editable', 'off', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.84 0.95 0.84]);
    ui.persistenceStatusField.Layout.Row = 1;
    ui.persistenceStatusField.Layout.Column = [1 2];
    ui.persistenceFailureField = uieditfield(persistenceStatusGrid, 'text', ...
        'Value', 'Failed to load the specified configuration', ...
        'Editable', 'off', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1.00 0.88 0.88], ...
        'Visible', 'off');
    ui.persistenceFailureField.Layout.Row = 1;
    ui.persistenceFailureField.Layout.Column = 1;
    ui.dismissPersistenceFailureButton = uibutton(persistenceStatusGrid, 'push', ...
        'Text', 'X', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.90 0.48 0.48], ...
        'Visible', 'off', ...
        'ButtonPushedFcn', @dismissPersistenceFailure);
    ui.dismissPersistenceFailureButton.Layout.Row = 1;
    ui.dismissPersistenceFailureButton.Layout.Column = 2;

    %% Lattice editing plane and view controls

    editorPanel = uipanel(mainGrid, ...
        'BorderType', 'none', ...
        'BackgroundColor', [1 1 1], ...
        'AutoResizeChildren', 'off');
    editorPanel.Layout.Row = 1;
    editorPanel.Layout.Column = 2;

    ui.axes = uiaxes(editorPanel, ...
        'Tag', 'LatticeAxes', ...
        'Box', 'on', ...
        'Color', [1 1 1], ...
        'XColor', [0.38 0.41 0.46], ...
        'YColor', [0.38 0.41 0.46], ...
        'GridColor', [0.82 0.84 0.87], ...
        'MinorGridColor', [0.90 0.91 0.93], ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'ButtonDownFcn', @axesClicked);
    xlabel(ui.axes, 'X', 'FontSize', axisLabelFontSize, 'Rotation', 0);
    ylabel(ui.axes, 'Y', 'FontSize', axisLabelFontSize, 'Rotation', 0);
    title(ui.axes, 'Lattice editing plane', ...
        'FontSize', axisTitleFontSize, 'FontWeight', 'normal');
    % Keep the visible plotting frame square without allowing MATLAB to
    % return either limit mode to automatic. In contrast, axis equal can
    % make the horizontal limits automatic when the UIAxes is rectangular.
    ui.axes.PlotBoxAspectRatio = [1 1 1];
    ui.axes.PlotBoxAspectRatioMode = 'manual';
    ui.axes.XLim = [-5 5];
    ui.axes.YLim = [-5 5];
    ui.axes.XLimMode = 'manual';
    ui.axes.YLimMode = 'manual';
    disableDefaultInteractivity(ui.axes);
    ui.axes.Toolbar.Visible = 'off';

    editorIconDirectory = fullfile(matlabroot, 'toolbox', 'matlab', ...
        'uitools', 'uifigureappjs', 'release', 'toolbox', 'matlab', ...
        'rich_text_component', 'src', 'js', 'rtc', 'images', ...
        'figurefloatingpalette');
    panDefaultIcon = fullfile(editorIconDirectory, 'pan_default.svg');
    panActiveIcon = fullfile(editorIconDirectory, 'pan_active.svg');
    homeIcon = fullfile(editorIconDirectory, 'restoreview_action.svg');
    ui.panButton = uibutton(editorPanel, 'state', ...
        'Text', '✋', ...
        'Tag', 'PanLatticeViewButton', ...
        'Tooltip', 'Pan editing plane', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'ValueChangedFcn', @panModeChanged);
    if isfile(panDefaultIcon)
        ui.panButton.Icon = panDefaultIcon;
        ui.panButton.Text = '';
    end
    ui.recenterButton = uibutton(editorPanel, 'push', ...
        'Text', '⌂', ...
        'Tag', 'RecenterLatticeViewButton', ...
        'FontWeight', 'bold', ...
        'Tooltip', 'Recenter editing plane', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'ButtonPushedFcn', @recenterEditingPlane);
    if isfile(homeIcon)
        ui.recenterButton.Icon = homeIcon;
        ui.recenterButton.Text = '';
    end

    % A single static web layer supplies true alpha compositing while the
    % blocking sublattice editor owns input. It is cheaper and more faithful
    % than capturing and repeatedly redrawing a dimmed image of this window.
    blockingOverlayHtml = [ ...
        '<!doctype html><html><head><meta charset="utf-8"><style>' ...
        'html,body{width:100%;height:100%;margin:0;overflow:hidden;' ...
        'background:rgba(160,165,172,0.45);}' ...
        'body{display:flex;align-items:center;justify-content:center;' ...
        'font:600 18px "Segoe UI",Arial,sans-serif;color:#454a52;' ...
        'user-select:none;}' ...
        '.message{padding:9px 22px;border:1px solid ' ...
        'rgba(110,115,122,0.45);border-radius:4px;' ...
        'background:rgba(250,250,250,0.88);}' ...
        '</style></head><body><div class="message">' ...
        'Sublattice editor is open</div></body></html>'];
    ui.sublatticeEditorOverlay = uihtml(mainGrid, ...
        'HTMLSource', blockingOverlayHtml, ...
        'Visible', 'off', ...
        'Tag', 'SublatticeEditorOverlay');
    ui.sublatticeEditorOverlay.Layout.Row = 1;
    ui.sublatticeEditorOverlay.Layout.Column = [1 2];

    %% Final layout and blocking editor loop

    applyTypography();
    editorPanel.SizeChangedFcn = @layoutEditor;
    drawnow;
    if ~isvalid(ui.figure)
        return;
    end
    layoutEditor();
    refreshAll();
    state.initialized = true;
    if state.cancelRequested
        delete(ui.figure);
        return;
    end
    ui.figure.Visible = 'on';
    drawnow;

    try
        uiwait(ui.figure);
    catch exception
        if ~isvalid(ui.figure)
            return;
        end
        rethrow(exception);
    end
    if isvalid(ui.figure)
        delete(ui.figure);
    end
    if ~configurationAccepted
        H_lattice = [];
        B = cell(0, 1);
        internalBlocks = cell(0, 1);
    elseif plotting
        launchModalViewer();
    end

    %% Typography and responsive layout

    function [panel,bodyHost,titleLabel,headerGrid] = ...
            createSectionPanel(parent, titleText, tag)
        % Build a fixed-height white header, separator, and flexible body.
        panel = uipanel(parent, ...
            'BackgroundColor', [1 1 1], ...
            'Tag', tag);
        sectionGrid = uigridlayout(panel, [3 1]);
        sectionGrid.RowHeight = {sectionHeaderHeight - 1, 1, '1x'};
        sectionGrid.ColumnWidth = {'1x'};
        sectionGrid.Padding = [0 0 0 0];
        sectionGrid.RowSpacing = 0;

        headerGrid = uigridlayout(sectionGrid, [1 2]);
        headerGrid.Layout.Row = 1;
        headerGrid.Layout.Column = 1;
        headerGrid.ColumnWidth = {'1x', '1x'};
        headerGrid.RowHeight = {'1x'};
        headerGrid.Padding = [7 0 7 0];
        headerGrid.ColumnSpacing = 5;
        headerGrid.BackgroundColor = [1 1 1];
        headerGrid.Tag = [tag 'Header'];
        titleLabel = uilabel(headerGrid, ...
            'Text', titleText, ...
            'FontSize', panelTitleFontSize, ...
            'VerticalAlignment', 'center');
        titleLabel.Layout.Row = 1;
        titleLabel.Layout.Column = [1 2];

        separator = uipanel(sectionGrid, ...
            'BorderType', 'none', ...
            'BackgroundColor', sectionSeparatorColor, ...
            'Tag', [tag 'Separator']);
        separator.Layout.Row = 2;
        separator.Layout.Column = 1;

        bodyHost = uigridlayout(sectionGrid, [1 1]);
        bodyHost.Layout.Row = 3;
        bodyHost.Layout.Column = 1;
        bodyHost.RowHeight = {'1x'};
        bodyHost.ColumnWidth = {'1x'};
        bodyHost.Padding = [0 0 0 0];
    end

    function applyTypography()
        % Normalize control fonts after construction, then restore role sizes.
        fontHandles = findall(ui.figure, '-property', 'FontSize');
        for fontIndex = 1:numel(fontHandles)
            fontHandles(fontIndex).FontSize = regularFontSize;
        end
        latticeTitleLabel.FontSize = panelTitleFontSize;
        ui.sublatticeTitleLabel.FontSize = panelTitleFontSize;
        controlsTitleLabel.FontSize = panelTitleFontSize;
        ui.formulaLabel.FontSize = regularFontSize + 1;
        ui.addButton.FontSize = buttonFontSize;
        ui.cancelButton.FontSize = buttonFontSize;
        ui.acceptButton.FontSize = buttonFontSize;
        ui.loadButton.FontSize = buttonFontSize;
        ui.saveButton.FontSize = buttonFontSize;
        ui.confirmDirectoryButton.FontSize = buttonFontSize;
        ui.cancelDirectoryButton.FontSize = buttonFontSize;
        ui.axes.FontSize = regularFontSize;
        ui.axes.XLabel.FontSize = axisLabelFontSize;
        ui.axes.YLabel.FontSize = axisLabelFontSize;
        ui.axes.Title.FontSize = axisTitleFontSize;
        ui.axes.Title.FontWeight = 'normal';
    end

    function layoutEditor(~, ~)
        % Fill the right column and keep both view buttons tied to the frame.
        if ~isvalid(editorPanel) || ~isvalid(ui.axes)
            return;
        end
        panelPosition = editorPanel.Position;
        width = max(panelPosition(3), 20);
        height = max(panelPosition(4), 24);
        ui.axes.Position = [10 10 max(width - 20, 1) max(height - 24, 1)];
        % UIAxes recomputes InnerPosition asynchronously after Position is
        % assigned. Flush that geometry before anchoring the view buttons so
        % they do not trail the plotting frame during an interactive resize.
        drawnow nocallbacks;
        if ~isvalid(editorPanel) || ~isvalid(ui.axes)
            return;
        end
        buttonSize = 26;
        buttonGap = 2;
        buttonFrameGap = 5;
        plotBox = axesPlotBoxPosition();
        plotRight = sum(plotBox([1 3]));
        plotTop = sum(plotBox([2 4]));
        pairWidth = 2*buttonSize + buttonGap;
        pairLeft = max(10, plotRight - pairWidth);
        panX = pairLeft;
        homeX = pairLeft + buttonSize + buttonGap;
        buttonY = max(10, min(height - buttonSize, plotTop + buttonFrameGap));
        ui.panButton.Position = [panX buttonY buttonSize buttonSize];
        ui.recenterButton.Position = [homeX buttonY buttonSize buttonSize];
    end

    function position = axesPlotBoxPosition()
        % Recover the square data frame inside UIAxes decorations and padding.
        position = ui.axes.InnerPosition;
        xSpan = diff(ui.axes.XLim);
        ySpan = diff(ui.axes.YLim);
        if xSpan <= 0 || ySpan <= 0
            return;
        end
        dataRatio = xSpan/ySpan;
        pixelRatio = position(3)/max(position(4), 1);
        if pixelRatio > dataRatio
            plotWidth = position(4)*dataRatio;
            position(1) = position(1) + (position(3) - plotWidth)/2;
            position(3) = plotWidth;
        elseif pixelRatio < dataRatio
            plotHeight = position(3)/dataRatio;
            position(2) = position(2) + (position(4) - plotHeight)/2;
            position(4) = plotHeight;
        end
    end

    %% Sublattice creation, editing, deletion, and properties

    function addSublattice(~, ~)
        % Accept a completed child model as a new rigid body at the origin.
        resetInteraction();
        [newH,newConfiguration] = openSublatticeEditor([]);
        if isempty(newConfiguration)
            return;
        end
        if isempty(newConfiguration.sites.uid)
            uialert(ui.figure, 'A sublattice must contain at least one site.', ...
                'Unable to add sublattice');
            return;
        end
        recordUndoState();
        newUid = state.nextUid;
        state.nextUid = state.nextUid + 1;
        state.sublattices(end + 1) = struct( ...
            'uid', newUid, ...
            'origin', [0 0], ...
            'rotationDegrees', 0, ...
            'configuration', newConfiguration, ...
            'H', newH);
        state.selectedUid = newUid;
        markConfigurationModified();
        refreshAll();
    end

    function editSublattice(~, ~)
        % Replace only the selected sublattice's internal model and H block.
        index = selectedSublatticeIndex();
        if isempty(index)
            return;
        end
        resetInteraction();
        [newH,newConfiguration] = openSublatticeEditor( ...
            state.sublattices(index).configuration);
        if isempty(newConfiguration)
            return;
        end
        if isempty(newConfiguration.sites.uid)
            uialert(ui.figure, 'A sublattice must contain at least one site.', ...
                'Unable to edit sublattice');
            return;
        end
        recordUndoState();
        state.sublattices(index).configuration = newConfiguration;
        state.sublattices(index).H = newH;
        markConfigurationModified();
        refreshAll();
    end

    function [newH,newConfiguration] = openSublatticeEditor(seed)
        % Keep the parent visible but inert while the blocking child runs.
        newH = [];
        newConfiguration = [];
        if ~isvalid(ui.figure)
            return;
        end
        outerFigure = ui.figure;
        blockingOverlay = ui.sublatticeEditorOverlay;
        overlayCleanup = onCleanup( ...
            @() hideSublatticeEditorOverlay(outerFigure, blockingOverlay));
        try
            showSublatticeEditorOverlay();
            if isempty(seed)
                [newH,newConfiguration] = editSublatticeConfiguration();
            else
                [newH,newConfiguration] = editSublatticeConfiguration(seed);
            end
        catch exception
            clear overlayCleanup;
            if isvalid(ui.figure)
                uialert(ui.figure, exception.message, ...
                    'Unable to open sublattice editor');
            end
            newH = [];
            newConfiguration = [];
            return;
        end
        clear overlayCleanup;
    end

    function showSublatticeEditorOverlay()
        % The overlay communicates modality and intercepts all parent input.
        if isvalid(ui.figure) && isvalid(ui.sublatticeEditorOverlay)
            ui.sublatticeEditorOverlay.Visible = 'on';
            drawnow;
        end
    end

    function deleteSelectedSublattice(~, ~)
        % Delete one complete rigid body and make the operation undoable.
        index = selectedSublatticeIndex();
        if isempty(index)
            return;
        end
        resetInteraction();
        recordUndoState();
        state.sublattices(index) = [];
        state.selectedUid = [];
        markConfigurationModified();
        refreshAll();
    end

    function originFieldChanged(source, ~, component)
        % Apply one Cartesian component without disturbing internal geometry.
        index = selectedSublatticeIndex();
        if isempty(index)
            return;
        end
        proposedOrigin = state.sublattices(index).origin;
        proposedOrigin(component) = source.Value;
        if any(~isfinite(proposedOrigin))
            source.Value = state.sublattices(index).origin(component);
            return;
        end
        if isequal(proposedOrigin, state.sublattices(index).origin)
            return;
        end
        recordUndoState();
        state.sublattices(index).origin = proposedOrigin;
        markConfigurationModified();
        refreshAll();
    end

    function rotationFieldChanged(source, ~)
        % Store an unwrapped angle; transforms use cosd/sind at render time.
        index = selectedSublatticeIndex();
        if isempty(index)
            return;
        end
        proposedRotation = source.Value;
        if ~isfinite(proposedRotation)
            source.Value = state.sublattices(index).rotationDegrees;
            return;
        end
        if proposedRotation == state.sublattices(index).rotationDegrees
            return;
        end
        recordUndoState();
        state.sublattices(index).rotationDegrees = proposedRotation;
        markConfigurationModified();
        refreshAll();
    end

    function index = selectedSublatticeIndex()
        % Resolve the stable UID and discard a selection deleted by undo/load.
        index = [];
        if isempty(state.selectedUid) || isempty(state.sublattices)
            return;
        end
        index = find([state.sublattices.uid] == state.selectedUid, 1, 'first');
        if isempty(index)
            state.selectedUid = [];
        end
    end

    function refreshSublatticeProperties()
        % Toggle the selection editor and publish the selected rigid transform.
        index = selectedSublatticeIndex();
        hasSelection = ~isempty(index);
        for controlIndex = 1:numel(ui.sublatticeSelectionControls)
            ui.sublatticeSelectionControls{controlIndex}.Visible = ...
                matlab.lang.OnOffSwitchState(hasSelection);
        end
        ui.sublatticeHint.Visible = ...
            matlab.lang.OnOffSwitchState(~hasSelection);
        if hasSelection
            ui.sublatticeTitleLabel.Text = ...
                sprintf('Sublattice %d properties', index);
            ui.xField.Value = state.sublattices(index).origin(1);
            ui.yField.Value = state.sublattices(index).origin(2);
            ui.rotationField.Value = state.sublattices(index).rotationDegrees;
            ui.xField.Enable = 'on';
            ui.yField.Enable = 'on';
            ui.rotationField.Enable = 'on';
            ui.editButton.Enable = 'on';
            ui.deleteButton.Enable = 'on';
        else
            ui.sublatticeTitleLabel.Text = 'Sublattice properties';
            ui.xField.Value = 0;
            ui.yField.Value = 0;
            ui.rotationField.Value = 0;
            ui.xField.Enable = 'off';
            ui.yField.Enable = 'off';
            ui.rotationField.Enable = 'off';
            ui.editButton.Enable = 'off';
            ui.deleteButton.Enable = 'off';
        end
    end

    %% Editing-plane interactions and navigation

    function axesClicked(~, event)
        % A primary click on empty space clears the current sublattice selection.
        if ui.panButton.Value
            beginPan(event);
            return;
        end
        if mouseButton(event) ~= 1
            return;
        end
        resetInteraction();
        state.selectedUid = [];
        refreshSublatticeProperties();
        redrawGraph();
    end

    function sublatticeClicked(~, event, uid)
        % Select by stable UID and stage one undo snapshot for a possible drag.
        if ui.panButton.Value
            beginPan(event);
            return;
        end
        if mouseButton(event) ~= 1
            return;
        end
        index = find([state.sublattices.uid] == uid, 1, 'first');
        if isempty(index)
            return;
        end
        state.selectedUid = uid;
        state.interaction = defaultInteractionState();
        dragUndoSnapshot = captureUndoState();
        dragUndoRecorded = false;
        state.interaction.mode = "drag_pending";
        state.interaction.startPoint = axesPoint(event);
        state.interaction.originalOrigin = state.sublattices(index).origin;
        state.interaction.originalCutoffDistance = state.cutoffDistance;
        state.interaction.uid = uid;
        refreshSublatticeProperties();
        redrawGraph();
    end

    function beginPan(event)
        % Pan in figure pixels so motion remains smooth outside the axes frame.
        if mouseButton(event) ~= 1
            return;
        end
        resetInteraction();
        state.interaction.mode = "pan_pending";
        state.interaction.startPoint = double(ui.figure.CurrentPoint(1, 1:2));
        state.interaction.originalXLimits = ui.axes.XLim;
        state.interaction.originalYLimits = ui.axes.YLim;
        ui.figure.Pointer = 'hand';
    end

    function pointerMoved(~, ~)
        % Advance the pending pan/drag state machine without recording zoom.
        mode = state.interaction.mode;
        if mode == "none"
            return;
        end
        if mode == "pan_pending" || mode == "pan"
            currentPixel = double(ui.figure.CurrentPoint(1, 1:2));
            pixelMovement = currentPixel - state.interaction.startPoint;
            if ~state.interaction.moved && hypot(pixelMovement(1), pixelMovement(2)) < 1
                return;
            end
            state.interaction.mode = "pan";
            state.interaction.moved = true;
            plotBox = axesPlotBoxPosition();
            xShift = pixelMovement(1)*diff(state.interaction.originalXLimits) ...
                /max(plotBox(3), 1);
            yShift = pixelMovement(2)*diff(state.interaction.originalYLimits) ...
                /max(plotBox(4), 1);
            ui.axes.XLim = state.interaction.originalXLimits - xShift;
            ui.axes.YLim = state.interaction.originalYLimits - yShift;
            ui.axes.XLimMode = 'manual';
            ui.axes.YLimMode = 'manual';
            drawnow limitrate nocallbacks;
            return;
        end
        if mode ~= "drag_pending" && mode ~= "drag"
            return;
        end
        index = find([state.sublattices.uid] == state.interaction.uid, 1, 'first');
        if isempty(index)
            resetInteraction();
            return;
        end
        currentPoint = axesPoint([]);
        movement = currentPoint - state.interaction.startPoint;
        threshold = 0.004*max(diff(ui.axes.XLim), diff(ui.axes.YLim));
        if ~state.interaction.moved && norm(movement) < threshold
            return;
        end
        state.interaction.mode = "drag";
        state.interaction.moved = true;
        proposedOrigin = state.interaction.originalOrigin + movement;
        if isequaln(proposedOrigin, state.sublattices(index).origin)
            return;
        end
        if ~dragUndoRecorded
            if isempty(dragUndoSnapshot)
                dragUndoSnapshot = captureUndoState();
            end
            pushUndoSnapshot(dragUndoSnapshot);
            dragUndoRecorded = true;
        end
        state.sublattices(index).origin = proposedOrigin;
        markConfigurationModified();
        refreshSublatticeProperties();
        % Use the cutoff at mouse-down as the baseline throughout a drag.
        % This lets the automatic floor follow d_min in either direction
        % without ratcheting upward through intermediate cursor positions.
        refreshGlobalMinimumDistance( ...
            state.interaction.originalCutoffDistance);
        redrawGraph();
        drawnow limitrate nocallbacks;
    end

    function pointerReleased(~, ~)
        % Pointer-up commits the already-applied motion by ending interaction.
        resetInteraction();
    end

    function resetInteraction()
        % Clear transient drag/pan data; semantic model state is untouched.
        state.interaction = defaultInteractionState();
        dragUndoSnapshot = [];
        dragUndoRecorded = false;
        if isfield(ui, 'figure') && isvalid(ui.figure)
            ui.figure.Pointer = 'arrow';
        end
    end

    function interaction = defaultInteractionState()
        % Keep all transient pointer fields in one predictable value object.
        interaction = struct( ...
            'mode', "none", ...
            'startPoint', [nan nan], ...
            'originalOrigin', [nan nan], ...
            'originalCutoffDistance', nan, ...
            'originalXLimits', [nan nan], ...
            'originalYLimits', [nan nan], ...
            'moved', false, ...
            'uid', []);
    end

    function panModeChanged(source, ~)
        % Reset any active gesture before changing the toggle's visual state.
        resetInteraction();
        if source.Value && isfile(panActiveIcon)
            source.Icon = panActiveIcon;
        elseif ~source.Value && isfile(panDefaultIcon)
            source.Icon = panDefaultIcon;
        end
    end

    function zoomEditingPlane(~, event)
        % Zoom about the pointer, bounded to practical data spans [0.5,200].
        if ~isvalid(ui.axes) || state.interaction.mode ~= "none"
            return;
        end
        hitObject = hittest(ui.figure);
        if isempty(hitObject)
            return;
        elseif isequal(hitObject, ui.axes)
            hitAxes = ui.axes;
        else
            hitAxes = ancestor(hitObject, 'axes');
        end
        if isempty(hitAxes) || ~isequal(hitAxes, ui.axes)
            return;
        end
        anchor = double(ui.axes.CurrentPoint(1, 1:2));
        oldXLimits = ui.axes.XLim;
        oldYLimits = ui.axes.YLim;
        if any(~isfinite(anchor)) ...
                || anchor(1) < oldXLimits(1) || anchor(1) > oldXLimits(2) ...
                || anchor(2) < oldYLimits(1) || anchor(2) > oldYLimits(2)
            return;
        end
        zoomSteps = double(event.VerticalScrollCount);
        if ~isfinite(zoomSteps) || zoomSteps == 0
            return;
        end
        currentSpan = diff(oldXLimits);
        requestedScale = (16/15)^zoomSteps;
        targetSpan = min(200, max(0.5, currentSpan*requestedScale));
        scale = targetSpan/currentSpan;
        ui.axes.XLim = anchor(1) + (oldXLimits - anchor(1))*scale;
        ui.axes.YLim = anchor(2) + (oldYLimits - anchor(2))*scale;
        ui.axes.XLimMode = 'manual';
        ui.axes.YLimMode = 'manual';
        drawnow limitrate nocallbacks;
    end

    function recenterEditingPlane(~, ~)
        % Restore the geometry-derived home view without entering undo history.
        resetInteraction();
        [xLimits,yLimits] = homeAxesLimits();
        ui.axes.XLim = xLimits;
        ui.axes.YLim = yLimits;
        ui.axes.XLimMode = 'manual';
        ui.axes.YLimMode = 'manual';
        drawnow limitrate nocallbacks;
    end

    function [xLimits,yLimits] = homeAxesLimits()
        % Fit all transformed sites into one square view with stable padding.
        positions = allGlobalPositions();
        if isempty(positions)
            xLimits = [-5 5];
            yLimits = [-5 5];
            return;
        end
        lower = min(positions, [], 1);
        upper = max(positions, [], 1);
        center = (lower + upper)/2;
        halfWidth = max(5, 0.60*max(upper - lower));
        xLimits = center(1) + [-halfWidth halfWidth];
        yLimits = center(2) + [-halfWidth halfWidth];
    end

    function handleKeyPress(sourceFigure, event)
        % Route undo/delete shortcuts while protecting active text entry.
        key = string(event.Key);
        modifiers = string(event.Modifier);
        isUndo = strcmpi(key, "z") ...
            && any(strcmpi(modifiers, "control") | strcmpi(modifiers, "command")) ...
            && ~any(strcmpi(modifiers, "shift"));
        if isUndo
            undoLastConfigurationChange();
            return;
        end
        isDelete = any(strcmpi(key, ["delete", "backspace"]));
        if isDelete && isempty(modifiers) && isequal(sourceFigure, ui.figure) ...
                && state.interaction.mode == "none" ...
                && ~keyboardFocusIsEditable(sourceFigure)
            deleteSelectedSublattice([], []);
        end
    end

    function tf = keyboardFocusIsEditable(figureHandle)
        % Numeric and text fields expose Editable; buttons and axes do not.
        focusedObject = figureHandle.CurrentObject;
        tf = ~isempty(focusedObject) && isvalid(focusedObject) ...
            && isprop(focusedObject, 'Editable') ...
            && strcmpi(string(focusedObject.Editable), "on");
    end

    function recordUndoState()
        % Capture state immediately before one semantic mutation.
        pushUndoSnapshot(captureUndoState());
    end

    function pushUndoSnapshot(snapshot)
        % Bound memory while retaining the most recent semantic edits.
        if isempty(snapshot)
            return;
        end
        undoStack{end + 1, 1} = snapshot;
        if numel(undoStack) > undoLimit
            undoStack(1) = [];
        end
    end

    function snapshot = captureUndoState()
        % Camera, selection, and persistence UI deliberately stay out of undo.
        snapshot = struct( ...
            'sublattices', state.sublattices, ...
            'nextUid', state.nextUid, ...
            'couplingProfile', state.couplingProfile, ...
            'gamma', state.gamma, ...
            'cutoffDistance', state.cutoffDistance);
    end

    function undoLastConfigurationChange()
        % Restore one model snapshot, retaining a still-valid current selection.
        if isempty(undoStack)
            return;
        end
        snapshot = undoStack{end};
        undoStack(end) = [];
        resetInteraction();
        state.sublattices = snapshot.sublattices;
        state.nextUid = snapshot.nextUid;
        state.couplingProfile = snapshot.couplingProfile;
        state.gamma = snapshot.gamma;
        state.cutoffDistance = snapshot.cutoffDistance;
        if ~isempty(state.selectedUid) ...
                && ~any([state.sublattices.uid] == state.selectedUid)
            state.selectedUid = [];
        end
        markConfigurationModified();
        refreshAll();
    end

    function point = axesPoint(event)
        % Prefer event intersection data; fall back to the axes current point.
        point = double(ui.axes.CurrentPoint(1, 1:2));
        if ~isempty(event) && isprop(event, 'IntersectionPoint')
            candidate = event.IntersectionPoint;
            if numel(candidate) >= 2 && all(isfinite(candidate(1:2)))
                point = double(candidate(1:2));
            end
        end
    end

    function button = mouseButton(event)
        % Normalize UI event and figure selection types to MATLAB button codes.
        button = 1;
        if ~isempty(event) && isprop(event, 'Button') && ~isempty(event.Button)
            button = event.Button;
        elseif strcmp(ui.figure.SelectionType, 'alt')
            button = 3;
        end
    end

    %% Coupling controls, model assembly, and geometry

    function profileChanged(source, ~)
        % Change only the automatic cross-sublattice distance law.
        newProfile = string(source.Value);
        if newProfile == state.couplingProfile
            return;
        end
        recordUndoState();
        state.couplingProfile = newProfile;
        markConfigurationModified();
        refreshAll();
    end

    function cutoffChanged(source, ~)
        % Enforce the global d_min floor before recording a cutoff edit.
        proposedCutoff = source.Value;
        if isnan(proposedCutoff) || proposedCutoff < 0
            source.Value = state.cutoffDistance;
            return;
        end
        proposedCutoff = enforceCutoffFloor( ...
            proposedCutoff, state.globalDMin);
        source.Value = proposedCutoff;
        if proposedCutoff == state.cutoffDistance
            return;
        end
        recordUndoState();
        state.cutoffDistance = proposedCutoff;
        markConfigurationModified();
        refreshAll();
    end

    function gammaChanged(source, ~)
        % The decay exponent/rate must remain finite and strictly positive.
        if ~isfinite(source.Value) || source.Value <= 0
            source.Value = state.gamma;
            return;
        end
        if source.Value == state.gamma
            return;
        end
        recordUndoState();
        state.gamma = source.Value;
        markConfigurationModified();
        refreshAll();
    end

    function formula = couplingFormula()
        % Keep the displayed LaTeX convention synchronized with evaluation.
        switch state.couplingProfile
            case "Power law"
                formula = ['$\kappa_{i,j}=(1+r_{i,j})^{-\gamma}' ...
                    '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
            case "Exponential"
                formula = ['$\kappa_{i,j}=\exp(-\gamma r_{i,j})' ...
                    '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
            case "Fiber asymptotic"
                formula = ['$\kappa_{i,j}=(1+r_{i,j})^{-1/2}\exp(-\gamma r_{i,j})' ...
                    '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
        end
    end

    function label = cutoffLabelText()
        % Show the current global normalization scale when it is defined.
        distance = state.globalDMin;
        if isfinite(distance)
            label = sprintf('Range cutoff (dₘᵢₙ=%.3f)', distance);
        else
            label = 'Range cutoff';
        end
    end

    function refreshGlobalMinimumDistance(cutoffBaseline)
        % Recompute transformed geometry and update the absolute cutoff floor.
        if nargin < 1
            cutoffBaseline = state.cutoffDistance;
        end
        state.globalDMin = globalMinimumDistance();
        state.cutoffDistance = enforceCutoffFloor( ...
            cutoffBaseline, state.globalDMin);
        ui.cutoffLabel.Text = cutoffLabelText();
        ui.cutoffField.Value = state.cutoffDistance;
    end

    function cutoff = enforceCutoffFloor(cutoff, distance)
        % Add one floating-point increment so a corrected cutoff exceeds d_min.
        if isfinite(distance) && distance > 0 && cutoff < distance
            cutoff = distance + eps(max(1, distance));
        end
    end

    function distance = globalMinimumDistance()
        % Measure every distinct pair after sublattice transforms are applied.
        distance = nan;
        positions = allGlobalPositions();
        if size(positions, 1) < 2
            return;
        end
        distance = Inf;
        for firstIndex = 1:size(positions, 1) - 1
            separations = vecnorm( ...
                positions(firstIndex + 1:end, :) - positions(firstIndex, :), ...
                2, 2);
            if ~isempty(separations)
                distance = min(distance, min(separations));
            end
        end
        if ~isfinite(distance)
            distance = nan;
        end
    end

    function positions = globalPositionsForSublattice(index)
        % Apply the rigid transform using row-vector coordinate storage.
        localPositions = state.sublattices(index).configuration.sites.position;
        angle = state.sublattices(index).rotationDegrees;
        rotation = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
        positions = localPositions*rotation.' + state.sublattices(index).origin;
    end

    function positions = allGlobalPositions()
        % Concatenation order is identical to H_lattice block ordering.
        positions = zeros(0, 2);
        for index = 1:numel(state.sublattices)
            positions = [positions; globalPositionsForSublattice(index)]; %#ok<AGROW>
        end
    end

    function [output,indices,blocks] = calculateCompositeHamiltonian()
        % Preserve internal blocks exactly, then fill reciprocal cross blocks.
        blocks = reshape({state.sublattices.H}, [], 1);
        indices = cell(numel(blocks), 1);
        if isempty(blocks)
            output = [];
            return;
        end
        counts = zeros(numel(blocks), 1);
        for blockIndex = 1:numel(blocks)
            block = blocks{blockIndex};
            count = numel(state.sublattices(blockIndex).configuration.sites.uid);
            if ~isnumeric(block) || ~isequal(size(block), [count count]) ...
                    || any(~isfinite(real(block)), 'all') ...
                    || any(~isfinite(imag(block)), 'all')
                error('lattice_config:InvalidSublatticeHamiltonian', ...
                    'Sublattice %d has an invalid internal Hamiltonian.', blockIndex);
            end
            counts(blockIndex) = count;
        end
        offsets = [0; cumsum(counts)];
        for blockIndex = 1:numel(blocks)
            indices{blockIndex} = offsets(blockIndex) + (1:counts(blockIndex));
        end
        output = blkdiag(blocks{:});
        if numel(blocks) < 2
            return;
        end

        dMin = globalMinimumDistance();
        if ~isfinite(dMin) || dMin < 1e-8
            error('lattice_config:CoincidentSites', ...
                ['Two lattice sites coincide. Edit or move a sublattice before ' ...
                 'building the composite Hamiltonian.']);
        end
        for firstIndex = 1:numel(blocks) - 1
            firstPositions = globalPositionsForSublattice(firstIndex);
            for secondIndex = firstIndex + 1:numel(blocks)
                secondPositions = globalPositionsForSublattice(secondIndex);
                distances = hypot( ...
                    secondPositions(:, 1) - firstPositions(:, 1).', ...
                    secondPositions(:, 2) - firstPositions(:, 2).');
                crossBlock = crossCouplingValues(distances, dMin);
                firstRange = offsets(firstIndex) + (1:counts(firstIndex));
                secondRange = offsets(secondIndex) + (1:counts(secondIndex));
                output(secondRange, firstRange) = crossBlock;
                output(firstRange, secondRange) = crossBlock.';
            end
        end
        if any(~isfinite(real(output)), 'all') || any(~isfinite(imag(output)), 'all')
            error('lattice_config:NonfiniteHamiltonian', ...
                'The composite Hamiltonian contains a nonfinite value.');
        end
    end

    function values = crossCouplingValues(distances, dMin)
        % Evaluate the selected profile only inside the absolute cutoff radius.
        values = zeros(size(distances));
        if isempty(distances) || ~isfinite(dMin) || dMin <= 0
            return;
        end
        included = distances <= state.cutoffDistance*(1 + 1e-10);
        q = distances(included)/dMin;
        q(abs(q - 1) <= 1e-10) = 1;
        switch state.couplingProfile
            case "Power law"
                values(included) = q.^(-state.gamma);
            case "Exponential"
                values(included) = exp(-state.gamma*(q - 1));
            case "Fiber asymptotic"
                values(included) = q.^(-1/2).*exp(-state.gamma*(q - 1));
        end
    end

    %% Graph rendering

    function refreshAll()
        % Synchronize derived geometry, controls, persistence state, and graph.
        refreshGlobalMinimumDistance();
        ui.profileDropdown.Value = char(state.couplingProfile);
        ui.formulaLabel.Text = couplingFormula();
        ui.cutoffField.Value = state.cutoffDistance;
        ui.gammaField.Value = state.gamma;
        refreshSublatticeProperties();
        refreshPersistenceControls();
        redrawGraph();
    end

    function redrawGraph()
        % Replace only tagged graph objects; axes configuration stays intact.
        if ~isvalid(ui.axes)
            return;
        end
        delete(findall(ui.axes, 'Tag', 'CompositeGraphObject'));
        hold(ui.axes, 'on');
        dMin = state.globalDMin;
        drawCrossSublatticeEdges(dMin);
        for sublatticeIndex = 1:numel(state.sublattices)
            drawOneSublattice(sublatticeIndex);
        end
        if isempty(state.sublattices)
            text(ui.axes, 0, 0, 'Add a sublattice to begin', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', regularFontSize + 4, ...
                'Color', [0.56 0.59 0.64], ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'Tag', 'CompositeGraphObject');
        end
        hold(ui.axes, 'off');
    end

    function drawCrossSublatticeEdges(dMin)
        % Batch straight reciprocal cross links by perceptual strength.
        if numel(state.sublattices) < 2 || ~isfinite(dMin) || dMin < 1e-8
            return;
        end
        strengths = [];
        segments = zeros(0, 4);
        for firstIndex = 1:numel(state.sublattices) - 1
            firstPositions = globalPositionsForSublattice(firstIndex);
            for secondIndex = firstIndex + 1:numel(state.sublattices)
                secondPositions = globalPositionsForSublattice(secondIndex);
                distances = hypot( ...
                    secondPositions(:, 1) - firstPositions(:, 1).', ...
                    secondPositions(:, 2) - firstPositions(:, 2).');
                values = crossCouplingValues(distances, dMin);
                [secondSites,firstSites] = find(values > 1e-14);
                for edgeIndex = 1:numel(firstSites)
                    firstPoint = firstPositions(firstSites(edgeIndex), :);
                    secondPoint = secondPositions(secondSites(edgeIndex), :);
                    segments(end + 1, :) = [firstPoint secondPoint]; %#ok<AGROW>
                    strengths(end + 1, 1) = values( ...
                        secondSites(edgeIndex), firstSites(edgeIndex)); %#ok<AGROW>
                end
            end
        end
        if isempty(segments)
            return;
        end
        % Quantized batches keep large cross-coupled arrays responsive. A
        % square-root display scale expands contrast among weak links, while
        % the near-white endpoint lets vanishing couplings become genuinely
        % faint instead of collecting at a conspicuous gray floor.
        normalizedStrengths = min(1, max(0, strengths));
        displayStrengths = sqrt(normalizedStrengths);
        binCount = 20;
        bins = min(binCount, floor(displayStrengths*binCount) + 1);
        for bin = 1:binCount
            edgeIndices = find(bins == bin);
            if isempty(edgeIndices)
                continue;
            end
            strength = mean(normalizedStrengths(edgeIndices));
            displayStrength = mean(displayStrengths(edgeIndices));
            shade = 0.96 - 0.74*displayStrength;
            xData = reshape([segments(edgeIndices, 1).'; ...
                segments(edgeIndices, 3).'; nan(1, numel(edgeIndices))], 1, []);
            yData = reshape([segments(edgeIndices, 2).'; ...
                segments(edgeIndices, 4).'; nan(1, numel(edgeIndices))], 1, []);
            line(ui.axes, xData, yData, ...
                'Color', [shade shade shade], ...
                'LineWidth', 0.6 + 1.2*strength, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'Tag', 'CompositeGraphObject');
        end
    end

    function drawOneSublattice(index)
        % Render a rigid sublattice with one stable UID-derived display color.
        sublattice = state.sublattices(index);
        positions = globalPositionsForSublattice(index);
        count = size(positions, 1);
        if count == 0
            return;
        end
        isSelected = isequal(state.selectedUid, sublattice.uid);
        uniformColor = sublatticeDisplayColor(sublattice.uid);
        drawInternalSublatticeEdges(sublattice, positions, uniformColor);

        lower = min(positions, [], 1);
        upper = max(positions, [], 1);
        span = max(upper - lower);
        padding = max(0.18, 0.06*max(span, 1));
        if isSelected
            rectangle(ui.axes, ...
                'Position', [lower - padding, upper - lower + 2*padding], ...
                'Curvature', 0.04, ...
                'EdgeColor', uniformColor, ...
                'LineStyle', ':', ...
                'LineWidth', 1.5, ...
                'FaceColor', 'none', ...
                'Tag', 'CompositeGraphObject', ...
                'ButtonDownFcn', @(source,event) ...
                    sublatticeClicked(source,event,sublattice.uid));
        end

        if isSelected
            scatter(ui.axes, positions(:, 1), positions(:, 2), 310, ...
                'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', uniformColor, ...
                'LineWidth', 1.8, ...
                'HitTest', 'off', ...
                'PickableParts', 'none', ...
                'Tag', 'CompositeGraphObject');
        end
        scatter(ui.axes, positions(:, 1), positions(:, 2), 190, ...
            'MarkerFaceColor', uniformColor, ...
            'MarkerEdgeColor', uniformColor, ...
            'LineWidth', 1, ...
            'Tag', 'CompositeGraphObject', ...
            'ButtonDownFcn', @(source,event) ...
                sublatticeClicked(source,event,sublattice.uid));
    end

    function color = sublatticeDisplayColor(uid)
        % A sublattice uses one hue for every internal node, line, and arrow.
        palette = [ ...
            0.12 0.47 0.71; ...
            0.85 0.33 0.10; ...
            0.17 0.63 0.17; ...
            0.58 0.40 0.74; ...
            0.55 0.34 0.29; ...
            0.89 0.47 0.76; ...
            0.38 0.38 0.38];
        color = palette(mod(uid - 1, size(palette, 1)) + 1, :);
    end

    function drawInternalSublatticeEdges(sublattice, positions, color)
        % Retain internal curves, dots, and asymmetry in a uniform outer color.
        configuration = sublattice.configuration;
        internalH = sublattice.H;
        count = size(internalH, 1);
        localDMin = localMinimumDistance(configuration.sites.position);
        manualMask = configuration.coupling.manualMask;
        ratios = configuration.coupling.syntheticDistanceRatio;
        for firstSite = 1:count - 1
            for secondSite = firstSite + 1:count
                forward = internalH(secondSite, firstSite);
                backward = internalH(firstSite, secondSite);
                forwardManual = manualMask(firstSite, secondSite);
                backwardManual = manualMask(secondSite, firstSite);
                if abs(forward) <= 1e-14 && abs(backward) <= 1e-14 ...
                        && ~forwardManual && ~backwardManual
                    continue;
                end
                comparisonScale = max(abs([forward backward]));
                asymmetric = comparisonScale > 0 ...
                    && abs(forward - backward) > 1e-10*comparisonScale;
                if ~asymmetric
                    directionalRatios = [ratios(firstSite, secondSite), ...
                        ratios(secondSite, firstSite)];
                    directionalRatios = directionalRatios(~isnan(directionalRatios));
                    if isempty(directionalRatios)
                        effectiveDistance = nan;
                    else
                        effectiveDistance = mean(directionalRatios)*localDMin;
                    end
                    drawInternalDirection(positions(firstSite, :), ...
                        positions(secondSite, :), effectiveDistance, 0, ...
                        false, color, sublattice.uid, max(abs([forward backward])));
                else
                    actualDistance = norm(positions(secondSite, :) ...
                        - positions(firstSite, :));
                    pathOffset = min(0.08, 0.07*actualDistance);
                    forwardDistance = manualEffectiveDistance( ...
                        ratios(firstSite, secondSite), localDMin);
                    backwardDistance = manualEffectiveDistance( ...
                        ratios(secondSite, firstSite), localDMin);
                    drawInternalDirection(positions(firstSite, :), ...
                        positions(secondSite, :), forwardDistance, pathOffset, ...
                        true, color, sublattice.uid, abs(forward));
                    drawInternalDirection(positions(secondSite, :), ...
                        positions(firstSite, :), backwardDistance, pathOffset, ...
                        true, color, sublattice.uid, abs(backward));
                end
            end
        end
    end

    function distance = localMinimumDistance(positions)
        % Internal synthetic ratios remain normalized to the local d_min.
        distance = nan;
        for firstSite = 1:size(positions, 1) - 1
            separations = vecnorm(positions(firstSite + 1:end, :) ...
                - positions(firstSite, :), 2, 2);
            pairMinimum = min(separations, [], 'omitnan');
            if ~isempty(pairMinimum) ...
                    && (isnan(distance) || pairMinimum < distance)
                distance = pairMinimum;
            end
        end
    end

    function distance = manualEffectiveDistance(ratio, localDMin)
        % Preserve NaN as the sentinel for unmodified straight geometry.
        if isnan(ratio) || isnan(localDMin)
            distance = nan;
        else
            distance = ratio*localDMin;
        end
    end

    function drawInternalDirection(source, destination, effectiveDistance, ...
            offset, showArrow, color, uid, magnitude)
        % Encode effective distance by curve length or a dotted shortcut.
        delta = destination - source;
        actualDistance = norm(delta);
        if actualDistance <= 1e-12
            return;
        end
        direction = delta/actualDistance;
        normal = [-direction(2) direction(1)];
        source = source + offset*normal;
        destination = destination + offset*normal;
        path = [source; destination];
        if isinf(effectiveDistance) && effectiveDistance > 0
            path = simplifiedCurvedPath(source, destination, 4*actualDistance);
        elseif isfinite(effectiveDistance) ...
                && effectiveDistance > actualDistance + 1e-6*max(1, actualDistance)
            path = simplifiedCurvedPath(source, destination, effectiveDistance);
        end

        if magnitude <= 1e-14
            drawInternalPath(path, ':', 1.0, color, uid);
        elseif isfinite(effectiveDistance) ...
                && actualDistance - effectiveDistance > 0.02*actualDistance
            solidFraction = max(0, min(1, effectiveDistance/actualDistance));
            firstBreak = solidFraction/2;
            secondBreak = 1 - firstBreak;
            drawInternalPath(pathFraction(path, 0, firstBreak), ...
                '-', 1.4, color, uid);
            drawInternalPath(pathFraction(path, firstBreak, secondBreak), ...
                ':', 1.6, color, uid);
            drawInternalPath(pathFraction(path, secondBreak, 1), ...
                '-', 1.4, color, uid);
        else
            drawInternalPath(path, '-', 1.4, color, uid);
        end
        if showArrow
            drawInternalArrow(path, color, uid, actualDistance);
        end
    end

    function drawInternalPath(path, style, width, color, uid)
        % Every internal path remains selectable as part of its rigid body.
        if size(path, 1) < 2
            return;
        end
        line(ui.axes, path(:, 1), path(:, 2), ...
            'Color', color, ...
            'LineStyle', style, ...
            'LineWidth', width, ...
            'Tag', 'CompositeGraphObject', ...
            'ButtonDownFcn', @(source,event) ...
                sublatticeClicked(source,event,uid));
    end

    function path = simplifiedCurvedPath(source, destination, effectiveDistance)
        % The outer preview favors a lightweight quadratic approximation.
        delta = destination - source;
        chord = norm(delta);
        direction = delta/chord;
        normal = [-direction(2) direction(1)];
        excess = max(0, effectiveDistance - chord);
        height = min(1.5*chord, sqrt(max(excess*chord, 0)));
        control = (source + destination)/2 + height*normal;
        parameter = linspace(0, 1, 41).';
        path = (1 - parameter).^2.*source ...
            + 2*(1 - parameter).*parameter.*control ...
            + parameter.^2.*destination;
    end

    function segment = pathFraction(path, firstFraction, secondFraction)
        % Extract a solid/dotted interval by accumulated arc length.
        parameter = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
        if parameter(end) <= eps
            segment = path([1 end], :);
            return;
        end
        parameter = parameter/parameter(end);
        samples = unique([firstFraction; ...
            parameter(parameter > firstFraction & parameter < secondFraction); ...
            secondFraction]);
        segment = [interp1(parameter, path(:, 1), samples), ...
            interp1(parameter, path(:, 2), samples)];
    end

    function drawInternalArrow(path, color, uid, actualDistance)
        % Use a compact data-space arrow suited to the simplified preview.
        [center,tangent] = pathPointAndTangent(path, 0.62);
        normal = [-tangent(2) tangent(1)];
        arrowLength = min(0.075, max(0.025, 0.07*actualDistance));
        arrowWidth = 0.65*arrowLength;
        tip = center + arrowLength*tangent;
        base = center - 0.45*arrowLength*tangent;
        vertices = [tip; base + arrowWidth*normal; base - arrowWidth*normal];
        patch(ui.axes, ...
            'XData', vertices(:, 1), ...
            'YData', vertices(:, 2), ...
            'FaceColor', color, ...
            'EdgeColor', color, ...
            'Tag', 'CompositeGraphObject', ...
            'ButtonDownFcn', @(source,event) ...
                sublatticeClicked(source,event,uid));
    end

    function [point,tangent] = pathPointAndTangent(path, fraction)
        % Sample a unit tangent at an arc-length fraction of a polyline.
        cumulative = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
        if cumulative(end) <= eps
            point = path(1, :);
            tangent = [1 0];
            return;
        end
        target = fraction*cumulative(end);
        segmentIndex = find(cumulative(2:end) >= target, 1, 'first');
        segmentLength = cumulative(segmentIndex + 1) - cumulative(segmentIndex);
        localFraction = (target - cumulative(segmentIndex))/max(segmentLength, eps);
        point = path(segmentIndex, :) + localFraction ...
            *(path(segmentIndex + 1, :) - path(segmentIndex, :));
        tangent = path(segmentIndex + 1, :) - path(segmentIndex, :);
        tangent = tangent/max(norm(tangent), eps);
    end

    %% Persistence workflow

    function beginPersistenceAction(action)
        % Replace Load/Save with an explicit directory confirmation step.
        state.pendingPersistenceAction = action;
        state.persistenceMode = "confirm";
        refreshPersistenceControls();
    end

    function cancelPersistenceAction(~, ~)
        % Abandon confirmation without changing either path field.
        state.pendingPersistenceAction = "";
        state.persistenceMode = "path";
        refreshPersistenceControls();
    end

    function confirmPersistenceAction(~, ~)
        % Execute the staged action only after reading both current path fields.
        action = state.pendingPersistenceAction;
        try
            directory = string(strtrim(ui.folderPathField.Value));
            fileName = string(strtrim(ui.fileNameField.Value));
            if action == "save"
                saveLatticeConfiguration(directory, fileName);
                showPersistenceSuccess("Saved");
            elseif action == "load"
                loadedConfiguration = readLatticeConfiguration(directory, fileName);
                preLoadSnapshot = captureUndoState();
                installConfigurationModel(loadedConfiguration);
                pushUndoSnapshot(preLoadSnapshot);
                resetInteraction();
                showPersistenceSuccess("Loaded");
                refreshAll();
                recenterEditingPlane([], []);
            end
        catch exception
            showPersistenceFailure(action, exception.message);
        end
    end

    function showPersistenceSuccess(message)
        % Collapse the path rows into the centered green status field.
        state.pendingPersistenceAction = "";
        state.persistenceMode = "success";
        ui.persistenceStatusField.Value = char(message);
        ui.persistenceStatusField.BackgroundColor = [0.84 0.95 0.84];
        refreshPersistenceControls();
    end

    function showPersistenceFailure(action, message)
        % Keep the concise inline failure while exposing details in an alert.
        state.pendingPersistenceAction = "";
        state.persistenceMode = "failure";
        ui.persistenceFailureField.Value = sprintf( ...
            'Failed to %s the specified configuration', char(action));
        refreshPersistenceControls();
        if isvalid(ui.figure)
            uialert(ui.figure, message, 'Lattice configuration error');
        end
    end

    function dismissPersistenceFailure(~, ~)
        % Restore editable folder and file-name rows after a failed operation.
        state.persistenceMode = "path";
        refreshPersistenceControls();
    end

    function persistenceLocationChanged(~, ~)
        % The edit fields own their text; only presentation mode is state.
        if state.persistenceMode ~= "confirm"
            state.persistenceMode = "path";
        end
        refreshPersistenceControls();
    end

    function markConfigurationModified()
        % A semantic edit clears transient Saved/Loaded/Failed presentation.
        if state.persistenceMode == "success" || state.persistenceMode == "failure"
            state.persistenceMode = "path";
        end
    end

    function refreshPersistenceControls()
        % Map the finite-state persistence mode onto mutually exclusive rows.
        mode = state.persistenceMode;
        persistenceButtonGrid.Visible = onOff(mode ~= "confirm");
        confirmationButtonGrid.Visible = onOff(mode == "confirm");
        persistenceInputGrid.Visible = onOff(mode == "path" || mode == "confirm");
        persistenceStatusGrid.Visible = onOff(mode == "success" || mode == "failure");
        ui.persistenceStatusField.Visible = onOff(mode == "success");
        ui.persistenceFailureField.Visible = onOff(mode == "failure");
        ui.dismissPersistenceFailureButton.Visible = onOff(mode == "failure");
    end

    function value = onOff(condition)
        % Convert logical expressions to graphics visibility strings.
        if condition
            value = 'on';
        else
            value = 'off';
        end
    end

    %% Composite configuration file I/O

    function filePath = configurationFilePath(directory, fileName)
        % The GUI separates directory and base name; storage uses one path.
        filePath = fullfile(char(directory), normalizedMatFileName(fileName));
    end

    function filePath = resolveSavedLatticePath(fileInput)
        % Normalize the headless-load argument and require an existing MAT-file.
        validText = (ischar(fileInput) && isrow(fileInput)) ...
            || (isstring(fileInput) && isscalar(fileInput));
        if ~validText
            error('lattice_config:InvalidLatticeFile', ...
                'The saved lattice file must be a character vector or string scalar.');
        end
        fileText = strtrim(string(fileInput));
        if ismissing(fileText) || strlength(fileText) == 0
            error('lattice_config:InvalidLatticeFile', ...
                'Enter the name or path of a saved lattice file.');
        end
        filePath = char(fileText);
        [~,~,extension] = fileparts(filePath);
        if isempty(extension)
            filePath = [filePath '.mat'];
        elseif ~strcmpi(extension, '.mat')
            error('lattice_config:InvalidLatticeFile', ...
                'A saved lattice file must use the .mat extension.');
        end
        if ~isfile(filePath)
            error('lattice_config:MissingLatticeFile', ...
                'Saved lattice file not found: %s', filePath);
        end
    end

    function normalizedName = normalizedMatFileName(fileName)
        % Reject paths, non-MAT extensions, and Windows-reserved base names.
        fileName = string(strtrim(fileName));
        if ~isscalar(fileName) || strlength(fileName) == 0
            error('lattice_config:InvalidFileName', 'Enter a file name.');
        end
        characters = char(fileName);
        if any(characters < 32) || any(ismember(characters, '<>:"/\|?*'))
            error('lattice_config:InvalidFileName', ...
                'The file name contains a character that Windows does not allow.');
        end
        [folderPart,baseName,extension] = fileparts(characters);
        if ~isempty(folderPart) || isempty(baseName) ...
                || (~isempty(extension) && ~strcmpi(extension, '.mat'))
            error('lattice_config:InvalidFileName', ...
                'Enter a file name with no extension, or use the .mat extension.');
        end
        reservedNames = ["CON","PRN","AUX","NUL", ...
            "COM1","COM2","COM3","COM4","COM5","COM6","COM7","COM8","COM9", ...
            "LPT1","LPT2","LPT3","LPT4","LPT5","LPT6","LPT7","LPT8","LPT9"];
        if any(strcmpi(string(baseName), reservedNames))
            error('lattice_config:InvalidFileName', ...
                'The file name is reserved by Windows.');
        end
        normalizedName = [baseName '.mat'];
    end

    function saveLatticeConfiguration(directory, fileName)
        % Save atomically through a same-directory temporary MAT-file.
        if ~isfolder(directory)
            error('lattice_config:MissingSaveDirectory', ...
                'The specified directory does not exist.');
        end
        latticeConfiguration = latticeConfigurationFromState();
        targetPath = configurationFilePath(directory, fileName);
        temporaryPath = [tempname(char(directory)) '.mat'];
        temporaryCleanup = onCleanup(@() deleteTemporaryConfiguration(temporaryPath));
        save(temporaryPath, 'latticeConfiguration', '-v7');
        [moveSucceeded,moveMessage] = movefile(temporaryPath, targetPath, 'f');
        if ~moveSucceeded
            error('lattice_config:SaveMoveFailed', ...
                'The lattice file could not be placed in the specified directory: %s', ...
                moveMessage);
        end
        clear temporaryCleanup;
    end

    function deleteTemporaryConfiguration(filePath)
        % onCleanup removes an incomplete temporary save after any exception.
        if isfile(filePath)
            delete(filePath);
        end
    end

    function configuration = latticeConfigurationFromState()
        % Serialize only semantic model data; UI, selection, and undo are absent.
        configuration = struct();
        % This token names the MAT schema rather than the public function.
        configuration.format = 'lattice_configuration';
        configuration.nextUid = state.nextUid;
        configuration.coupling = struct( ...
            'profile', char(state.couplingProfile), ...
            'gamma', state.gamma, ...
            'cutoffDistance', state.cutoffDistance);
        configuration.sublattices = state.sublattices;
    end

    function configuration = readLatticeConfiguration(directory, fileName)
        % Adapt the two GUI fields to the common path-based reader.
        configuration = readLatticeConfigurationFile( ...
            configurationFilePath(directory, fileName));
    end

    function configuration = readLatticeConfigurationFile(filePath)
        % Load only the expected variable, then validate before installation.
        if ~isfile(filePath)
            error('lattice_config:MissingLatticeFile', ...
                'Saved lattice file not found: %s', filePath);
        end
        loadedData = load(filePath, 'latticeConfiguration');
        if ~isfield(loadedData, 'latticeConfiguration')
            error('lattice_config:MissingSavedConfiguration', ...
                'The MAT-file does not contain a lattice configuration.');
        end
        configuration = validateLatticeConfiguration(loadedData.latticeConfiguration);
    end

    function configuration = validateLatticeConfiguration(configuration)
        % Enforce the exact outer schema and canonical numeric shapes/types.
        requireConfiguration(isscalar(configuration) && isstruct(configuration), ...
            'The saved lattice configuration must be a scalar structure.');
        requireFields(configuration, {'format','nextUid','coupling','sublattices'});
        requireConfiguration(strcmp(string(configuration.format), ...
            "lattice_configuration"), ...
            'The MAT-file is not a lattice configuration.');

        coupling = configuration.coupling;
        requireConfiguration(isscalar(coupling) && isstruct(coupling), ...
            'The saved coupling properties are invalid.');
        requireFields(coupling, {'profile','gamma','cutoffDistance'});
        profile = string(coupling.profile);
        requireConfiguration(isscalar(profile) && any(profile == ...
            ["Power law","Exponential","Fiber asymptotic"]), ...
            'The saved coupling profile is invalid.');
        gamma = coupling.gamma;
        requireConfiguration(isnumeric(gamma) && isreal(gamma) ...
            && isscalar(gamma) && isfinite(gamma) && gamma > 0, ...
            'The saved decay parameter is invalid.');
        cutoff = coupling.cutoffDistance;
        requireConfiguration(isnumeric(cutoff) && isreal(cutoff) ...
            && isscalar(cutoff) && ~isnan(cutoff) && cutoff >= 0, ...
            'The saved range cutoff is invalid.');

        sublattices = configuration.sublattices;
        requireConfiguration(isstruct(sublattices), ...
            'The saved sublattice array is invalid.');
        if ~isempty(sublattices)
            requireFields(sublattices, ...
                {'uid','origin','rotationDegrees','configuration','H'});
        end
        uids = zeros(numel(sublattices), 1);
        canonical = struct('uid', {}, 'origin', {}, 'rotationDegrees', {}, ...
            'configuration', {}, 'H', {});
        for index = 1:numel(sublattices)
            item = sublattices(index);
            requireConfiguration(isnumeric(item.uid) && isreal(item.uid) ...
                && isscalar(item.uid) && isfinite(item.uid) ...
                && item.uid >= 1 && item.uid == fix(item.uid), ...
                'A saved sublattice identifier is invalid.');
            requireConfiguration(isnumeric(item.origin) && isreal(item.origin) ...
                && numel(item.origin) == 2 && all(isfinite(item.origin)), ...
                'A saved sublattice origin is invalid.');
            requireConfiguration(isnumeric(item.rotationDegrees) ...
                && isreal(item.rotationDegrees) && isscalar(item.rotationDegrees) ...
                && isfinite(item.rotationDegrees), ...
                'A saved sublattice rotation is invalid.');
            validateEmbeddedSublatticeConfiguration(item.configuration, item.H);
            uids(index) = item.uid;
            canonical(index) = struct( ...
                'uid', double(item.uid), ...
                'origin', reshape(double(item.origin), 1, 2), ...
                'rotationDegrees', double(item.rotationDegrees), ...
                'configuration', item.configuration, ...
                'H', double(item.H));
        end
        requireConfiguration(numel(unique(uids)) == numel(uids), ...
            'Saved sublattice identifiers must be unique.');
        nextUid = configuration.nextUid;
        minimumNextUid = 1;
        if ~isempty(uids)
            minimumNextUid = max(uids) + 1;
        end
        requireConfiguration(isnumeric(nextUid) && isreal(nextUid) ...
            && isscalar(nextUid) && isfinite(nextUid) ...
            && nextUid == fix(nextUid) && nextUid >= minimumNextUid, ...
            'The saved next-sublattice identifier is invalid.');

        configuration.nextUid = double(nextUid);
        configuration.coupling = struct( ...
            'profile', char(profile), ...
            'gamma', double(gamma), ...
            'cutoffDistance', double(cutoff));
        configuration.sublattices = canonical;
    end

    function validateEmbeddedSublatticeConfiguration(configuration, internalH)
        % Check the fields required to place and preserve each internal block.
        requireConfiguration(isscalar(configuration) && isstruct(configuration), ...
            'An embedded sublattice configuration is invalid.');
        requireFields(configuration, ...
            {'format','nextUid','sites','coupling','disorder','preset'});
        requireConfiguration(strcmp(string(configuration.format), ...
            "sublattice_configuration"), ...
            'An embedded record is not a sublattice configuration.');
        requireConfiguration(isscalar(configuration.sites) ...
            && isstruct(configuration.sites), ...
            'Embedded sublattice site data are invalid.');
        requireFields(configuration.sites, {'uid','position','detuning'});
        uid = configuration.sites.uid;
        positions = configuration.sites.position;
        detuning = configuration.sites.detuning;
        count = numel(uid);
        requireConfiguration(count > 0 && isnumeric(uid) && isreal(uid) ...
            && iscolumn(uid) && all(isfinite(uid)) ...
            && all(uid >= 1) && all(uid == fix(uid)) ...
            && numel(unique(uid)) == count, ...
            'Embedded sublattice site identifiers are invalid.');
        requireConfiguration(isnumeric(positions) && isreal(positions) ...
            && isequal(size(positions), [count 2]) ...
            && all(isfinite(positions), 'all'), ...
            'Embedded sublattice coordinates are invalid.');
        requireConfiguration(isnumeric(detuning) && isreal(detuning) ...
            && iscolumn(detuning) && numel(detuning) == count ...
            && all(isfinite(detuning)), ...
            'Embedded sublattice detunings are invalid.');
        requireConfiguration(isnumeric(internalH) ...
            && isequal(size(internalH), [count count]) ...
            && all(isfinite(real(internalH)), 'all') ...
            && all(isfinite(imag(internalH)), 'all'), ...
            'An embedded sublattice Hamiltonian is invalid.');
    end

    function requireFields(value, names)
        % Exact field sets prevent silent acceptance of misspelled properties.
        requireConfiguration(numel(fieldnames(value)) == numel(names), ...
            'The saved lattice configuration contains unexpected fields.');
        for nameIndex = 1:numel(names)
            requireConfiguration(isfield(value, names{nameIndex}), ...
                sprintf('The saved lattice configuration is missing "%s".', ...
                names{nameIndex}));
        end
    end

    function requireConfiguration(condition, message)
        % Give every schema failure one stable public error identifier.
        if ~condition
            error('lattice_config:InvalidLatticeConfiguration', ...
                '%s', message);
        end
    end

    function installConfigurationModel(configuration)
        % Install validated semantics, then derive d_min and the cutoff floor.
        state.nextUid = configuration.nextUid;
        state.sublattices = configuration.sublattices;
        state.selectedUid = [];
        state.couplingProfile = string(configuration.coupling.profile);
        state.gamma = configuration.coupling.gamma;
        state.cutoffDistance = configuration.coupling.cutoffDistance;
        state.globalDMin = globalMinimumDistance();
        state.cutoffDistance = enforceCutoffFloor( ...
            state.cutoffDistance, state.globalDMin);
    end

    %% Build, modal-viewer launch, and shutdown

    function launchModalViewer()
        % Plotting is ancillary: a viewer failure must not discard a valid H.
        try
            systems = modalViewerSystems( ...
                H_lattice, internalBlocks, state.sublattices);
            openModalViewer(systems);
        catch exception
            warning('lattice_config:ModalViewerFailed', ...
                'The Hamiltonian was built, but its modal viewer failed: %s', ...
                exception.message);
        end
    end

    function acceptConfiguration(~, ~)
        % Assemble once, then resume the blocking editor only on success.
        ui.acceptButton.Enable = 'off';
        drawnow;
        try
            [H_lattice,B,internalBlocks] = calculateCompositeHamiltonian();
        catch exception
            ui.acceptButton.Enable = 'on';
            uialert(ui.figure, exception.message, 'Unable to build lattice');
            return;
        end
        configurationAccepted = true;
        if isvalid(ui.figure) && strcmp(ui.figure.WaitStatus, 'waiting')
            uiresume(ui.figure);
        end
    end

    function cancelConfiguration(~, ~)
        % Cancel returns the documented empty outputs at every construction stage.
        H_lattice = [];
        B = cell(0, 1);
        internalBlocks = cell(0, 1);
        configurationAccepted = false;
        if ~state.initialized
            state.cancelRequested = true;
            return;
        end
        if isvalid(ui.figure) && strcmp(ui.figure.WaitStatus, 'waiting')
            uiresume(ui.figure);
        elseif isvalid(ui.figure)
            delete(ui.figure);
        end
    end
end

function systems = modalViewerSystems(H_lattice, internalBlocks, sublattices)
%MODALVIEWERSYSTEMS Freeze Hamiltonians and graph geometry for the viewer.
%   The modal viewer deliberately owns no editor handles. Local coordinates
%   are retained for an individual sublattice, while the composite entry
%   uses the same translated and rotated site order as H_LATTICE.

    systemTemplate = struct( ...
        'label', '', ...
        'H', complex(zeros(0)), ...
        'positions', zeros(0, 2), ...
        'primitives', modalViewerEmptyPrimitives());
    sublatticeCount = numel(sublattices);
    systems = repmat(systemTemplate, sublatticeCount + 1, 1);

    globalPositions = cell(sublatticeCount, 1);
    globalPrimitives = modalViewerEmptyPrimitives();
    siteOffsets = zeros(sublatticeCount + 1, 1);
    for sublatticeIndex = 1:sublatticeCount
        sublattice = sublattices(sublatticeIndex);
        localPositions = double(sublattice.configuration.sites.position);
        localPrimitives = modalViewerInternalPrimitives( ...
            sublattice.configuration, internalBlocks{sublatticeIndex});

        rotation = [ ...
            cosd(sublattice.rotationDegrees), -sind(sublattice.rotationDegrees); ...
            sind(sublattice.rotationDegrees),  cosd(sublattice.rotationDegrees)];
        globalPositions{sublatticeIndex} = localPositions*rotation.' ...
            + reshape(sublattice.origin, 1, 2);
        transformedPrimitives = modalViewerTransformPrimitives( ...
            localPrimitives, rotation, sublattice.origin);
        globalPrimitives = [globalPrimitives transformedPrimitives]; %#ok<AGROW>

        systems(sublatticeIndex + 1) = struct( ...
            'label', sprintf('Subsystem %d', sublatticeIndex), ...
            'H', internalBlocks{sublatticeIndex}, ...
            'positions', localPositions, ...
            'primitives', localPrimitives);
        siteOffsets(sublatticeIndex + 1) = siteOffsets(sublatticeIndex) ...
            + size(localPositions, 1);
    end

    % Inter-sublattice blocks are reciprocal and geometric, so one neutral
    % straight primitive represents both directions in modal-profile plots.
    crossEdgeCount = 0;
    for firstSublattice = 1:sublatticeCount - 1
        firstIndices = siteOffsets(firstSublattice) + 1:siteOffsets(firstSublattice + 1);
        for secondSublattice = firstSublattice + 1:sublatticeCount
            secondIndices = siteOffsets(secondSublattice) ...
                + 1:siteOffsets(secondSublattice + 1);
            forwardBlock = H_lattice(secondIndices, firstIndices);
            backwardBlock = H_lattice(firstIndices, secondIndices).';
            crossEdgeCount = crossEdgeCount ...
                + nnz(abs(forwardBlock) > 1e-14 | abs(backwardBlock) > 1e-14);
        end
    end
    crossPrimitives = repmat( ...
        modalViewerPrimitive(zeros(2), '-', 0.8), 1, crossEdgeCount);
    crossPrimitiveIndex = 0;
    for firstSublattice = 1:sublatticeCount - 1
        for secondSublattice = firstSublattice + 1:sublatticeCount
            firstPositions = globalPositions{firstSublattice};
            secondPositions = globalPositions{secondSublattice};
            for firstSite = 1:size(firstPositions, 1)
                firstGlobalIndex = siteOffsets(firstSublattice) + firstSite;
                for secondSite = 1:size(secondPositions, 1)
                    secondGlobalIndex = siteOffsets(secondSublattice) + secondSite;
                    forward = H_lattice(secondGlobalIndex, firstGlobalIndex);
                    backward = H_lattice(firstGlobalIndex, secondGlobalIndex);
                    if abs(forward) <= 1e-14 && abs(backward) <= 1e-14
                        continue;
                    end
                    crossPrimitiveIndex = crossPrimitiveIndex + 1;
                    crossPrimitives(crossPrimitiveIndex) = modalViewerPrimitive( ...
                        [firstPositions(firstSite, :); ...
                         secondPositions(secondSite, :)], '-', 0.8);
                end
            end
        end
    end
    globalPrimitives = [globalPrimitives crossPrimitives];

    if isempty(globalPositions)
        compositePositions = zeros(0, 2);
    else
        compositePositions = vertcat(globalPositions{:});
    end
    systems(1) = struct( ...
        'label', 'Full system', ...
        'H', H_lattice, ...
        'positions', compositePositions, ...
        'primitives', globalPrimitives);
end

function openModalViewer(systems)
%OPENMODALVIEWER Open a maximized, nonblocking modal-analysis window.
%   Eigensystems are solved on first selection and then cached. Switching
%   between spectrum, intensity, and phase therefore redraws only graphics.

    viewer = uifigure( ...
        'Name', 'DNLSE modal viewer', ...
        'Tag', 'LatticeModalViewerFigure', ...
        'Color', [0.97 0.98 0.99], ...
        'Visible', 'off');
    viewer.WindowState = 'maximized';

    viewerGrid = uigridlayout(viewer, [2 1]);
    viewerGrid.RowHeight = {54, '1x'};
    viewerGrid.ColumnWidth = {'1x'};
    viewerGrid.Padding = [8 8 8 8];
    viewerGrid.RowSpacing = 7;

    toolbar = uipanel(viewerGrid, ...
        'BackgroundColor', [1 1 1], ...
        'BorderColor', [0.72 0.72 0.72]);
    toolbar.Layout.Row = 1;
    toolbar.Layout.Column = 1;
    toolbarGrid = uigridlayout(toolbar, [1 5]);
    toolbarGrid.ColumnWidth = {104, 260, 56, 220, '1x'};
    toolbarGrid.RowHeight = {'1x'};
    toolbarGrid.Padding = [12 7 12 7];
    toolbarGrid.ColumnSpacing = 8;

    hamiltonianLabel = uilabel(toolbarGrid, ...
        'Text', 'Hamiltonian', ...
        'FontSize', 16, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'center');
    hamiltonianLabel.Layout.Column = 1;
    systemLabels = {systems.label};
    hamiltonianDropdown = uidropdown(toolbarGrid, ...
        'Items', systemLabels, ...
        'ItemsData', 1:numel(systems), ...
        'Value', 1, ...
        'FontSize', 16, ...
        'Tag', 'ModalHamiltonianDropdown');
    hamiltonianDropdown.Layout.Column = 2;

    viewLabel = uilabel(toolbarGrid, ...
        'Text', 'View', ...
        'FontSize', 16, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'center');
    viewLabel.Layout.Column = 3;
    viewDropdown = uidropdown(toolbarGrid, ...
        'Items', {'Spectrum', 'Modal intensities', 'Modal phases'}, ...
        'Value', 'Spectrum', ...
        'FontSize', 16, ...
        'Tag', 'ModalViewDropdown');
    viewDropdown.Layout.Column = 4;

    plotPanel = uipanel(viewerGrid, ...
        'BackgroundColor', [1 1 1], ...
        'BorderType', 'none', ...
        'Tag', 'ModalPlotPanel');
    plotPanel.Layout.Row = 2;
    plotPanel.Layout.Column = 1;

    eigensystemCache = cell(numel(systems), 1);
    hamiltonianDropdown.ValueChangedFcn = @selectionChanged;
    viewDropdown.ValueChangedFcn = @selectionChanged;
    viewer.Visible = 'on';
    drawnow;
    redrawViewer();

    function selectionChanged(~, ~)
        % Either dropdown selection invalidates only the visible plot layout.
        redrawViewer();
    end

    function redrawViewer()
        % Lazily solve the selected H, then rebuild the requested plot tiles.
        if ~isvalid(viewer) || ~isvalid(plotPanel)
            return;
        end
        hamiltonianDropdown.Enable = 'off';
        viewDropdown.Enable = 'off';
        viewer.Pointer = 'watch';
        drawnow;
        cleanup = onCleanup(@() modalViewerRestoreControls( ...
            viewer, hamiltonianDropdown, viewDropdown));
        delete(plotPanel.Children);

        systemIndex = hamiltonianDropdown.Value;
        selectedSystem = systems(systemIndex);
        if isempty(selectedSystem.H)
            modalViewerMessage(plotPanel, 'No modes available');
            return;
        end
        if isempty(eigensystemCache{systemIndex})
            [psi, epsilon] = eig(selectedSystem.H, 'vector');
            [~, order] = sortrows([real(epsilon), imag(epsilon)], [1 2]);
            eigensystemCache{systemIndex} = struct( ...
                'epsilon', epsilon(order), ...
                'psi', psi(:, order));
        end
        solution = eigensystemCache{systemIndex};
        switch string(viewDropdown.Value)
            case "Spectrum"
                modalViewerSpectrum(plotPanel, solution.epsilon);
            case "Modal intensities"
                modalViewerProfiles(plotPanel, selectedSystem, ...
                    abs(solution.psi).^2, false);
            case "Modal phases"
                modalViewerProfiles(plotPanel, selectedSystem, ...
                    angle(solution.psi), true);
        end
        clear cleanup;
    end

end

function modalViewerRestoreControls(viewer, hamiltonianDropdown, viewDropdown)
%MODALVIEWERRESTORECONTROLS Restore controls after a redraw or error.
    if isvalid(viewer)
        viewer.Pointer = 'arrow';
    end
    if isvalid(hamiltonianDropdown)
        hamiltonianDropdown.Enable = 'on';
    end
    if isvalid(viewDropdown)
        viewDropdown.Enable = 'on';
    end
end

function modalViewerMessage(parent, message)
%MODALVIEWERMESSAGE Fill the plot panel with a centered status message.
    layout = uigridlayout(parent, [1 1]);
    layout.Padding = [0 0 0 0];
    label = uilabel(layout, ...
        'Text', message, ...
        'FontSize', 20, ...
        'FontColor', [0.45 0.48 0.52], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'center');
    label.Layout.Row = 1;
    label.Layout.Column = 1;
end

function modalViewerSpectrum(parent, epsilon)
%MODALVIEWERSPECTRUM Plot real and imaginary spectral parts side by side.
    modeCount = numel(epsilon);
    layout = tiledlayout(parent, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    realAxes = nexttile(layout);
    plot(realAxes, 1:modeCount, real(epsilon), 'k.', 'MarkerSize', 18)
    axis(realAxes, 'padded')
    xlim(realAxes, [0.5 modeCount + 0.5])
    title(realAxes, 'Real part', 'FontWeight', 'normal')
    ylabel(realAxes, 'Re(\epsilon)')
    box(realAxes, 'on')

    imaginaryAxes = nexttile(layout);
    plot(imaginaryAxes, 1:modeCount, imag(epsilon), 'k.', 'MarkerSize', 18)
    axis(imaginaryAxes, 'padded')
    xlim(imaginaryAxes, [0.5 modeCount + 0.5])
    title(imaginaryAxes, 'Imaginary part', 'FontWeight', 'normal')
    ylabel(imaginaryAxes, 'Im(\epsilon)')
    box(imaginaryAxes, 'on')

    title(layout, 'Eigenvalues', 'FontWeight', 'normal')
    xlabel(layout, 'Mode index')
end

function modalViewerProfiles(parent, system, values, phasePlot)
%MODALVIEWERPROFILES Draw every mode on the selected lattice geometry.
    modeCount = size(values, 2);
    columnCount = ceil(sqrt(modeCount));
    rowCount = ceil(modeCount/columnCount);
    layout = tiledlayout(parent, rowCount, columnCount, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    if phasePlot
        colorMap = ot_colormap('hsv',256);
        colorLimits = [-pi pi];
        colorbarLabel = 'Phase (rad)';
    else
        fullHotMap = ot_colormap('hot',256);
        colorMap = interp1(linspace(0, 1, 256), fullHotMap, ...
            linspace(0, 0.95, 256));
        % Values are abs(psi).^2. Stop at 95 percent of hot so the highest
        % intensity remains yellow instead of reaching the map's white tip.
        maximumIntensity = max(values, [], 'all');
        if maximumIntensity <= eps
            maximumIntensity = 1;
        end
        colorLimits = [0 maximumIntensity];
        colorbarLabel = 'Intensity';
    end

    [xLimits, yLimits] = modalViewerGraphLimits( ...
        system.positions, system.primitives);
    lastAxes = [];
    for modeIndex = 1:modeCount
        axesHandle = nexttile(layout);
        lastAxes = axesHandle;
        hold(axesHandle, 'on')
        for primitiveIndex = 1:numel(system.primitives)
            primitive = system.primitives(primitiveIndex);
            line(axesHandle, primitive.path(:, 1), primitive.path(:, 2), ...
                'Color', [0.68 0.68 0.68], ...
                'LineStyle', primitive.style, ...
                'LineWidth', primitive.width)
        end
        scatter(axesHandle, system.positions(:, 1), system.positions(:, 2), 72, ...
            values(:, modeIndex), 'filled', ...
            'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 0.5)
        colormap(axesHandle, colorMap)
        clim(axesHandle, colorLimits)
        xlim(axesHandle, xLimits)
        ylim(axesHandle, yLimits)
        axis(axesHandle, 'equal')
        axis(axesHandle, 'off')
        title(axesHandle, sprintf('Mode %d', modeIndex), ...
            'FontWeight', 'normal')
        hold(axesHandle, 'off')
    end

    colorbarHandle = colorbar(lastAxes);
    colorbarHandle.Layout.Tile = 'east';
    colorbarHandle.Label.String = colorbarLabel;
    if phasePlot
        colorbarHandle.Ticks = [-pi 0 pi];
        colorbarHandle.TickLabels = {'-\pi', '0', '\pi'};
    end
end

function primitives = modalViewerInternalPrimitives(configuration, H)
%MODALVIEWERINTERNALPRIMITIVES Preserve edited internal edge geometry.
%   H contains effective disorder values. The configuration retains the
%   underlying manual/automatic coupling state needed to reconstruct dots,
%   curves, and separated paths without restoring editor-specific colors.

    primitives = modalViewerEmptyPrimitives();
    positions = double(configuration.sites.position);
    count = size(positions, 1);
    if count < 2
        return;
    end
    dMin = modalViewerMinimumDistance(positions);
    baseCouplings = modalViewerBaseCouplings(configuration, positions, dMin);
    effectiveCouplings = H.';
    manualMask = logical(configuration.coupling.manualMask);

    for firstIndex = 1:count - 1
        for secondIndex = firstIndex + 1:count
            firstBaseValue = baseCouplings(firstIndex, secondIndex);
            secondBaseValue = baseCouplings(secondIndex, firstIndex);
            firstManual = manualMask(firstIndex, secondIndex);
            secondManual = manualMask(secondIndex, firstIndex);
            if abs(firstBaseValue) <= 1e-14 && abs(secondBaseValue) <= 1e-14 ...
                    && ~firstManual && ~secondManual
                continue;
            end

            firstValue = effectiveCouplings(firstIndex, secondIndex);
            secondValue = effectiveCouplings(secondIndex, firstIndex);
            comparisonScale = max(abs([firstValue secondValue]));
            valuesMatch = comparisonScale == 0 || ...
                abs(firstValue - secondValue) <= 1e-12*comparisonScale;
            if valuesMatch
                if string(configuration.disorder.type) == "Coupling disorder"
                    distanceRatio = mean([ ...
                        modalViewerDisplayedDistanceRatio(configuration, ...
                            baseCouplings, firstIndex, secondIndex, firstValue), ...
                        modalViewerDisplayedDistanceRatio(configuration, ...
                            baseCouplings, secondIndex, firstIndex, secondValue)]);
                else
                    distanceRatio = nan;
                    manualRatios = [ ...
                        configuration.coupling.syntheticDistanceRatio( ...
                            firstIndex, secondIndex), ...
                        configuration.coupling.syntheticDistanceRatio( ...
                            secondIndex, firstIndex)];
                    manualRatios = manualRatios(~isnan(manualRatios));
                    if ~isempty(manualRatios)
                        distanceRatio = mean(manualRatios);
                    end
                end
                newPrimitives = modalViewerDirectionalPrimitives( ...
                    positions(firstIndex, :), positions(secondIndex, :), ...
                    firstValue, modalViewerScaledDistance(distanceRatio, dMin), 0);
                primitives = [primitives newPrimitives]; %#ok<AGROW>
            else
                firstDistanceRatio = modalViewerDisplayedDistanceRatio( ...
                    configuration, baseCouplings, ...
                    firstIndex, secondIndex, firstValue);
                secondDistanceRatio = modalViewerDisplayedDistanceRatio( ...
                    configuration, baseCouplings, ...
                    secondIndex, firstIndex, secondValue);
                firstPrimitives = modalViewerDirectionalPrimitives( ...
                    positions(firstIndex, :), positions(secondIndex, :), ...
                    firstValue, modalViewerScaledDistance(firstDistanceRatio, dMin), ...
                    0.13);
                secondPrimitives = modalViewerDirectionalPrimitives( ...
                    positions(secondIndex, :), positions(firstIndex, :), ...
                    secondValue, modalViewerScaledDistance(secondDistanceRatio, dMin), ...
                    0.13);
                primitives = [primitives firstPrimitives secondPrimitives]; %#ok<AGROW>
            end
        end
    end
end

function baseCouplings = modalViewerBaseCouplings(configuration, positions, dMin)
%MODALVIEWERBASECOUPLINGS Rebuild pre-disorder source-to-destination values.
    count = size(positions, 1);
    baseCouplings = complex(zeros(count));
    if count >= 2 && isfinite(dMin)
        profile = string(configuration.coupling.profile);
        gamma = configuration.coupling.gamma;
        cutoffDistance = configuration.coupling.cutoffDistance;
        for firstIndex = 1:count - 1
            for secondIndex = firstIndex + 1:count
                separation = norm(positions(secondIndex, :) ...
                    - positions(firstIndex, :));
                if separation > cutoffDistance*(1 + 1e-10)
                    value = 0;
                else
                    distanceRatio = separation/dMin;
                    if abs(distanceRatio - 1) <= 1e-10
                        distanceRatio = 1;
                    end
                    value = modalViewerProfileValue( ...
                        distanceRatio, profile, gamma);
                end
                baseCouplings(firstIndex, secondIndex) = value;
                baseCouplings(secondIndex, firstIndex) = value;
            end
        end
    end
    manualMask = logical(configuration.coupling.manualMask);
    baseCouplings(manualMask) = configuration.coupling.manualValue(manualMask);
end

function value = modalViewerProfileValue(distanceRatio, profile, gamma)
%MODALVIEWERPROFILEVALUE Evaluate a normalized geometric coupling profile.
    switch string(profile)
        case "Power law"
            value = distanceRatio^(-gamma);
        case "Exponential"
            value = exp(-gamma*(distanceRatio - 1));
        case "Fiber asymptotic"
            value = distanceRatio^(-0.5)*exp(-gamma*(distanceRatio - 1));
        otherwise
            error('lattice_config:UnknownCouplingProfile', ...
                'Unknown coupling profile "%s".', profile);
    end
end

function ratio = modalViewerDisplayedDistanceRatio( ...
        configuration, baseCouplings, sourceIndex, destinationIndex, effectiveValue)
%MODALVIEWERDISPLAYEDDISTANCERATIO Include multiplicative disorder in paths.
    ratio = configuration.coupling.syntheticDistanceRatio( ...
        sourceIndex, destinationIndex);
    if string(configuration.disorder.type) ~= "Coupling disorder"
        return;
    end

    baseMagnitude = abs(baseCouplings(sourceIndex, destinationIndex));
    effectiveMagnitude = abs(effectiveValue);
    profile = string(configuration.coupling.profile);
    gamma = configuration.coupling.gamma;
    if isnan(ratio)
        ratio = modalViewerInferDistanceRatio(effectiveMagnitude, profile, gamma);
        return;
    end
    if baseMagnitude == effectiveMagnitude || baseMagnitude == 0
        return;
    end

    baseProfileRatio = modalViewerInferDistanceRatio( ...
        baseMagnitude, profile, gamma);
    effectiveProfileRatio = modalViewerInferDistanceRatio( ...
        effectiveMagnitude, profile, gamma);
    if isinf(effectiveProfileRatio)
        ratio = Inf;
    elseif isfinite(baseProfileRatio) && isfinite(ratio)
        ratio = max(0, ratio + effectiveProfileRatio - baseProfileRatio);
    end
end

function ratio = modalViewerInferDistanceRatio(magnitude, profile, gamma)
%MODALVIEWERINFERDISTANCERATIO Invert magnitude to q = d_eff/d_min.
    if magnitude == 0
        ratio = Inf;
        return;
    end
    switch string(profile)
        case "Power law"
            ratio = exp(-log(magnitude)/gamma);
        case "Exponential"
            ratio = max(0, 1 - log(magnitude)/gamma);
        case "Fiber asymptotic"
            ratio = modalViewerInverseFiberDistanceRatio(magnitude, gamma);
        otherwise
            error('lattice_config:UnknownCouplingProfile', ...
                'Unknown coupling profile "%s".', profile);
    end
end

function ratio = modalViewerInverseFiberDistanceRatio(magnitude, gamma)
%MODALVIEWERINVERSEFIBERDISTANCERATIO Invert the fiber profile in log(q).
    logMagnitude = log(magnitude);
    if logMagnitude == 0
        ratio = 1;
        return;
    elseif logMagnitude > 0
        lowerLogRatio = -2*(logMagnitude + 1);
        upperLogRatio = 0;
    else
        lowerLogRatio = 0;
        upperLogRatio = 2*(-logMagnitude + 1);
    end
    for iteration = 1:80
        middleLogRatio = (lowerLogRatio + upperLogRatio)/2;
        residual = gamma*expm1(middleLogRatio) ...
            + 0.5*middleLogRatio + logMagnitude;
        if residual < 0
            lowerLogRatio = middleLogRatio;
        else
            upperLogRatio = middleLogRatio;
        end
    end
    ratio = exp((lowerLogRatio + upperLogRatio)/2);
end

function distance = modalViewerScaledDistance(distanceRatio, dMin)
%MODALVIEWERSCALEDDISTANCE Convert normalized synthetic distance to length.
    if isnan(distanceRatio) || isnan(dMin)
        distance = nan;
    else
        distance = distanceRatio*dMin;
    end
end

function primitives = modalViewerDirectionalPrimitives( ...
        source, destination, value, syntheticDistance, offset)
%MODALVIEWERDIRECTIONALPRIMITIVES Build neutral solid/dotted path pieces.
    primitives = modalViewerEmptyPrimitives();
    delta = destination - source;
    actualDistance = norm(delta);
    if actualDistance <= 1e-12
        return;
    end
    direction = delta/actualDistance;
    normal = [-direction(2) direction(1)];
    source = source + offset*normal;
    destination = destination + offset*normal;
    path = [source; destination];

    if abs(value) <= 1e-14
        primitives = modalViewerPrimitive(path, ':', 1.0);
        return;
    end
    if isnan(syntheticDistance)
        effectiveDistance = actualDistance;
    elseif isinf(syntheticDistance)
        effectiveDistance = min(4*actualDistance, 24);
    else
        effectiveDistance = syntheticDistance;
    end

    tolerance = 1e-6*max(1, actualDistance);
    if effectiveDistance > actualDistance + tolerance
        path = modalViewerCurvedPath(source, destination, effectiveDistance);
    end
    isShortened = actualDistance - effectiveDistance > 0.02*actualDistance;
    if ~isShortened
        primitives = modalViewerPrimitive(path, '-', 0.8);
        return;
    end

    solidFraction = max(0, min(1, effectiveDistance/actualDistance));
    firstBreak = solidFraction/2;
    secondBreak = 1 - firstBreak;
    primitives(1) = modalViewerPrimitive( ...
        modalViewerPathFraction(path, 0, firstBreak), '-', 0.8);
    primitives(2) = modalViewerPrimitive( ...
        modalViewerPathFraction(path, firstBreak, secondBreak), ':', 1.0);
    primitives(3) = modalViewerPrimitive( ...
        modalViewerPathFraction(path, secondBreak, 1), '-', 0.8);
    primitives(4) = modalViewerPrimitive( ...
        modalViewerShortcutBoundary(path, firstBreak, actualDistance), '-', 1.0);
    primitives(5) = modalViewerPrimitive( ...
        modalViewerShortcutBoundary(path, secondBreak, actualDistance), '-', 1.0);
end

function primitive = modalViewerPrimitive(path, style, width)
%MODALVIEWERPRIMITIVE Package one neutral line primitive.
    primitive = struct('path', path, 'style', style, 'width', width);
end

function primitives = modalViewerEmptyPrimitives()
%MODALVIEWEREMPTYPRIMITIVES Return an empty typed primitive array.
    primitives = struct('path', {}, 'style', {}, 'width', {});
end

function output = modalViewerTransformPrimitives(primitives, rotation, origin)
%MODALVIEWERTRANSFORMPRIMITIVES Move local paths into composite coordinates.
    output = primitives;
    for primitiveIndex = 1:numel(output)
        output(primitiveIndex).path = output(primitiveIndex).path*rotation.' ...
            + reshape(origin, 1, 2);
    end
end

function boundary = modalViewerShortcutBoundary(path, fraction, actualDistance)
%MODALVIEWERSHORTCUTBOUNDARY Draw a small tick at a solid/dotted transition.
    [point, tangent] = modalViewerPathPointAndTangent(path, fraction);
    normal = [-tangent(2) tangent(1)];
    halfLength = 0.04*actualDistance;
    boundary = [point - halfLength*normal; point + halfLength*normal];
end

function path = modalViewerCurvedPath(source, destination, targetLength)
%MODALVIEWERCURVEDPATH Solve a quadratic arc with the requested length.
    chord = norm(destination - source);
    targetLength = max(chord, min(targetLength, max(4*chord, 24)));
    direction = (destination - source)/chord;
    normal = [-direction(2) direction(1)];
    low = 0;
    high = max(chord, 0.5);
    while modalViewerQuadraticArcLength(chord, high) < targetLength && high < 100
        high = 2*high;
    end
    for iteration = 1:45
        midpoint = (low + high)/2;
        if modalViewerQuadraticArcLength(chord, midpoint) < targetLength
            low = midpoint;
        else
            high = midpoint;
        end
    end
    sag = (low + high)/2;
    parameter = linspace(0, 1, 81).';
    straightPart = (1 - parameter).*source + parameter.*destination;
    path = straightPart + (4*sag*parameter.*(1 - parameter)).*normal;
end

function lengthValue = modalViewerQuadraticArcLength(chord, sag)
%MODALVIEWERQUADRATICARCLENGTH Return the exact quadratic-arc length.
    if sag <= eps*max(1, chord)
        lengthValue = chord;
        return;
    end
    lengthValue = 0.5*sqrt(chord^2 + 16*sag^2) ...
        + chord^2/(8*sag)*asinh(4*sag/chord);
end

function segment = modalViewerPathFraction(path, startFraction, endFraction)
%MODALVIEWERPATHFRACTION Extract an interval by path arc length.
    startFraction = max(0, min(1, startFraction));
    endFraction = max(startFraction, min(1, endFraction));
    cumulative = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
    total = cumulative(end);
    if total <= eps || endFraction - startFraction <= eps
        segment = zeros(0, 2);
        return;
    end
    startDistance = startFraction*total;
    endDistance = endFraction*total;
    interior = cumulative > startDistance & cumulative < endDistance;
    startPoint = modalViewerInterpolatePath(path, cumulative, startDistance);
    endPoint = modalViewerInterpolatePath(path, cumulative, endDistance);
    segment = [startPoint; path(interior, :); endPoint];
end

function point = modalViewerInterpolatePath(path, cumulative, distance)
%MODALVIEWERINTERPOLATEPATH Interpolate one point by path distance.
    if distance <= 0
        point = path(1, :);
        return;
    elseif distance >= cumulative(end)
        point = path(end, :);
        return;
    end
    upper = find(cumulative >= distance, 1, 'first');
    lower = upper - 1;
    fraction = (distance - cumulative(lower))/ ...
        max(cumulative(upper) - cumulative(lower), eps);
    point = path(lower, :) + fraction*(path(upper, :) - path(lower, :));
end

function [point, tangent] = modalViewerPathPointAndTangent(path, fraction)
%MODALVIEWERPATHPOINTANDTANGENT Sample a point and unit tangent by arc length.
    cumulative = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
    if cumulative(end) <= eps
        point = path(1, :);
        tangent = [1 0];
        return;
    end
    target = max(0, min(1, fraction))*cumulative(end);
    if target <= 0
        lower = 1;
        upper = 2;
    elseif target >= cumulative(end)
        lower = size(path, 1) - 1;
        upper = size(path, 1);
    else
        upper = find(cumulative >= target, 1, 'first');
        lower = upper - 1;
    end
    segmentLength = cumulative(upper) - cumulative(lower);
    localFraction = (target - cumulative(lower))/max(segmentLength, eps);
    point = path(lower, :) + localFraction*(path(upper, :) - path(lower, :));
    tangent = path(upper, :) - path(lower, :);
    tangent = tangent/max(norm(tangent), eps);
end

function dMin = modalViewerMinimumDistance(positions)
%MODALVIEWERMINIMUMDISTANCE Return the shortest nonzero site separation.
    if size(positions, 1) < 2
        dMin = nan;
        return;
    end
    dMin = Inf;
    for firstIndex = 1:size(positions, 1) - 1
        separations = vecnorm(positions(firstIndex + 1:end, :) ...
            - positions(firstIndex, :), 2, 2);
        separations = separations(separations > 1e-12);
        if ~isempty(separations)
            dMin = min(dMin, min(separations));
        end
    end
    if ~isfinite(dMin)
        dMin = nan;
    end
end

function [xLimits, yLimits] = modalViewerGraphLimits(positions, primitives)
%MODALVIEWERGRAPHLIMITS Center every mode in one shared square frame.
    graphPoints = positions;
    for primitiveIndex = 1:numel(primitives)
        graphPoints = [graphPoints; primitives(primitiveIndex).path]; %#ok<AGROW>
    end
    if isempty(graphPoints)
        xLimits = [-1 1];
        yLimits = [-1 1];
        return;
    end
    center = (min(graphPoints, [], 1) + max(graphPoints, [], 1))/2;
    coordinateSpan = max(max(graphPoints, [], 1) - min(graphPoints, [], 1));
    if coordinateSpan <= eps
        halfWidth = 1;
    else
        halfWidth = 0.58*coordinateSpan;
    end
    xLimits = center(1) + [-halfWidth halfWidth];
    yLimits = center(2) + [-halfWidth halfWidth];
end

function hideSublatticeEditorOverlay(figureHandle, overlayHandle)
%HIDESUBLATTICEEDITOROVERLAY Restore the parent after a blocking child editor.
    if isvalid(figureHandle) && isvalid(overlayHandle)
        overlayHandle.Visible = 'off';
        drawnow;
    end
end

function deleteFigureIfValid(figureHandle)
%DELETEFIGUREIFVALID Close an editor left behind by an interrupted workflow.
    if isvalid(figureHandle)
        delete(figureHandle);
    end
end

function html = justifiedControlsHtml(paragraph)
%JUSTIFIEDCONTROLSHTML Format static help copy with full justification.
    html = [ ...
        '<!doctype html><html><head><meta charset="utf-8"><style>' ...
        'html,body{width:100%;height:100%;margin:0;overflow:hidden;' ...
        'background:transparent;}' ...
        'p{margin:0;color:#262626;font:14px/1.25 "Segoe UI",Arial,sans-serif;' ...
        'text-align:justify;text-align-last:left;hyphens:none;}' ...
        '</style></head><body><p>' char(paragraph) '</p></body></html>'];
end

function position = defaultEditorWindowPosition( ...
        topLeft, leftColumnWidth, leftRowHeights, leftRowSpacing, ...
        mainPadding, mainColumnSpacing)
%DEFAULTEDITORWINDOWPOSITION Size a two-column editor from its left stack.
%   The plot column follows the total menu height, so adding vertical room
%   to the properties panel also enlarges the default square editing plane.

    leftStackHeight = sum(leftRowHeights) + ...
        max(0, numel(leftRowHeights) - 1)*leftRowSpacing;
    figureHeight = leftStackHeight + mainPadding(2) + mainPadding(4);

    % UIAxes spends slightly more vertical space on its title and X label
    % than horizontal space on its Y label and tick labels. This small
    % compensation places the rendered square frame at the width/height knee
    % without wasted column space; the result still scales with menu height.
    axesDecorationCompensation = 20;
    editorColumnWidth = max(300, leftStackHeight - axesDecorationCompensation);
    figureWidth = mainPadding(1) + leftColumnWidth + mainColumnSpacing + ...
        editorColumnWidth + mainPadding(3);

    position = [topLeft(1), topLeft(2) - figureHeight, ...
        figureWidth, figureHeight];
end

function [H,configuration] = editSublatticeConfiguration(seedConfiguration)
%EDITSUBLATTICECONFIGURATION Run the private editor for one sublattice.
%   This file-local editor supplies the Add and Edit workflows of
%   LATTICE_CONFIG. It returns the exact internal Hamiltonian H and
%   the complete editable CONFIGURATION record used by the outer editor.
%
%   SEEDCONFIGURATION may contain the current embedded record for an existing
%   sublattice. Omitting it creates a fresh interactive sublattice.
%
%   LATTICE DATA
%   Each site has a position (x,y), an on-site detuning, and a persistent
%   identifier. Display indices are assigned row by row: left to right,
%   then top to bottom. Internally, a directional coupling is stored as
%   coupling(source,destination); the returned matrix therefore satisfies
%
%       H(destination,source) = coupling(source,destination).
%
%   Site detunings occupy the diagonal of H. Reciprocal couplings need not
%   be equal, and each direction may carry an independent magnitude and
%   phase. Phases are entered in degrees and retained without wrapping.
%
%   GEOMETRIC COUPLINGS
%   Let r_ij = d_ij/d_min - 1, where d_min is the shortest nonzero current
%   site separation.
%   The selectable coupling profiles are
%
%       Power law:        kappa_ij = (1+r_ij)^(-gamma)
%       Exponential:      kappa_ij = exp(-gamma*r_ij)
%       Fiber asymptotic: kappa_ij = (1+r_ij)^(-1/2)*exp(-gamma*r_ij).
%
%   The absolute range cutoff removes automatic geometric edges at larger
%   coordinate separations; Inf disables truncation. A finite cutoff cannot
%   lie below the current d_min and is raised just above d_min when needed.
%   Otherwise, geometry edits preserve the chosen absolute radius. Editing
%   an edge creates a directional manual override containing a base complex
%   coupling and a profile-equivalent distance ratio. These values remain
%   fixed until Restore geometric couplings is pressed. The displayed
%   effective distance is measured in the same coordinate units as site
%   positions.
%
%   PRESET LATTICES
%   Square, triangular, honeycomb, two-dimensional J_x, Lieb, and ring
%   presets are included. Sublattice properties always expose the common
%   profile, cutoff, decay, and disorder controls. For a named preset with
%   no active element selection, the next panel shows family size, spacing,
%   and sublattice detunings. Selecting a site or edge temporarily exposes
%   its element fields; editing one first converts the generated lattice to
%   Custom while preserving the visible lattice and current camera view.
%
%   DISORDER
%   Disorder is an overlay on the editable base lattice. Anderson disorder
%   adds independent uniform samples in [-S/2,S/2] to the detunings.
%   Coupling disorder multiplies both directions of each unordered pair by
%   1+S*xi, with xi uniform in [-1/2,1/2]. Aubry–André disorder adds
%   S*cos(2*pi*beta*i), beta=(sqrt(5)-1)/2, using the displayed site index.
%   While a random disorder family remains active, existing samples stay
%   attached to site identifiers and unordered site pairs; new identifiers
%   receive new samples. Leaving and later re-entering a random family
%   starts a new realization. Property fields always show and modify the
%   underlying disorder-free values.
%
%   EDITOR CONTROLS
%   Double-click blank space to add a site. Left-click sites or edges to
%   select them; Shift-click and box selection select multiple sites. Drag
%   selected sites to move them, right-click a site to remove it, and use
%   Delete or Backspace to remove selected sites when no text field is
%   active. Ctrl+Z (or Command+Z) undoes sublattice edits. The hand button
%   enables panning, the mouse wheel zooms about the pointer, and the home
%   button restores the family view. Camera changes are not part of undo.
%
    %% Constants and editor state

    narginchk(0, 1);
    H = [];
    configuration = [];

    regularFontSize = 14;
    panelTitleFontSize = 22;
    buttonFontSize = 16;
    axisLabelFontSize = 14;
    axisTitleFontSize = 22;
    sectionHeaderHeight = 38;
    sectionSeparatorColor = [0.55 0.55 0.55];
    presetParameterFieldWidth = 72;
    readOnlyBackgroundColor = [0.94 0.94 0.94];
    readOnlyFontColor = [0.38 0.38 0.38];

    % State contains only the editable physical model and interaction data.
    % Disorder samples are stored separately and are applied on display/build.
    state = struct();
    state.nodes = struct('uid', {}, 'x', {}, 'y', {}, 'detuning', {});
    state.couplings = complex(zeros(0));       % Base value (source,destination).
    state.manual = false(0);                   % True for a directional override.
    state.manualPhase = nan(0);                % Unwrapped radians for the editor.
    state.syntheticDistanceRatio = nan(0);     % Manual profile-equivalent d_eff/d_min.
    state.nextUid = 1;
    state.gamma = 1;
    state.cutoffDistance = Inf;
    % Preset defaults follow preset spacing until the cutoff is edited directly.
    state.cutoffTracksPresetDefault = false;
    state.couplingProfile = "Exponential";
    state.disorder = defaultDisorderState();
    state.preset = "Custom";
    state.presetSettings = defaultPresetSettings();
    state.presetGlobals = defaultPresetGlobals();
    state.customSnapshot = [];
    state.selectedSiteUids = [];
    state.selectedEdgeUids = [];
    state.lastBlankClickTimer = [];
    state.lastBlankClickPosition = [nan nan];
    % Pending modes distinguish a click from a drag before any mutation.
    state.interaction = struct( ...
        'mode', "none", ...
        'startPoint', [nan nan], ...
        'currentPoint', [nan nan], ...
        'dragUids', [], ...
        'originalPositions', zeros(0, 2), ...
        'originalXLimits', [nan nan], ...
        'originalYLimits', [nan nan], ...
        'moved', false, ...
        'additive', false, ...
        'initialSelection', [], ...
        'clickedUid', []);
    state.initialized = false;
    state.cancelRequested = false;
    seededConfiguration = [];
    if nargin == 1
        seededConfiguration = validateSublatticeConfiguration(seedConfiguration);
    end

    % Undo holds semantic model snapshots only; view and selection stay live.
    undoStack = cell(0, 1);
    undoLimit = 40;
    % A drag records its pre-drag snapshot once, on the first accepted motion.
    dragUndoSnapshot = [];
    dragUndoRecorded = false;
    dragCutoffDistance = [];
    viewHomeHalfWidth = 5;

    %% Main window and column layout

    leftColumnWidth = 398;
    leftRowHeights = [191 157 355 44];
    leftRowSpacing = 5;
    mainPadding = [6 6 6 6];
    mainColumnSpacing = 6;
    defaultFigurePosition = defaultEditorWindowPosition( ...
        [617 866], leftColumnWidth, leftRowHeights, leftRowSpacing, ...
        mainPadding, mainColumnSpacing);

    ui = struct();
    ui.figure = uifigure( ...
        'Name', 'DNLSE sublattice configuration', ...
        'Tag', 'SublatticeConfigurationFigure', ...
        'Position', defaultFigurePosition, ...
        'WindowStyle', 'modal', ...
        'Color', [0.96 0.97 0.98], ...
        'Visible', 'off', ...
        'WindowButtonMotionFcn', @pointerMoved, ...
        'WindowButtonUpFcn', @pointerReleased, ...
        'WindowScrollWheelFcn', @zoomEditingPlane, ...
        'WindowKeyPressFcn', @handleUndoShortcut, ...
        'CloseRequestFcn', @cancelConfiguration);
    editorFigure = ui.figure;
    editorFigureCleanup = onCleanup(@() deleteFigureIfValid(editorFigure));

    mainGrid = uigridlayout(ui.figure, [1 2]);
    mainGrid.ColumnWidth = {leftColumnWidth, '1x'};
    mainGrid.RowHeight = {'1x'};
    mainGrid.Padding = mainPadding;
    mainGrid.ColumnSpacing = mainColumnSpacing;

    leftGrid = uigridlayout(mainGrid, [4 1]);
    leftGrid.Layout.Row = 1;
    leftGrid.Layout.Column = 1;
    leftGrid.Scrollable = 'off';
    leftGrid.RowHeight = num2cell(leftRowHeights);
    leftGrid.ColumnWidth = {'1x'};
    leftGrid.Padding = [0 0 0 0];
    leftGrid.RowSpacing = leftRowSpacing;

    %% Sublattice property controls

    % The preset selector shares the custom header with the coupling and
    % disorder controls, which remain visible for every preset.
    [sublatticePropertiesPanel,sublatticePropertiesBody, ...
        sublatticePropertiesLabel,propertyHeaderGrid] = ...
        createSectionPanel(leftGrid, 'Sublattice properties', ...
        'SublatticePropertiesPanel');
    sublatticePropertiesPanel.Layout.Row = 1;
    propertyHeaderGrid.ColumnWidth = {'1x', 114};
    sublatticePropertiesLabel.Tag = 'SublatticePropertiesHeaderLabel';
    sublatticePropertiesLabel.Layout.Column = 1;
    presetHeaderCell = uigridlayout(propertyHeaderGrid, [1 1]);
    presetHeaderCell.Layout.Row = 1;
    presetHeaderCell.Layout.Column = 2;
    presetHeaderCell.RowHeight = {'1x'};
    presetHeaderCell.ColumnWidth = {'1x'};
    presetHeaderCell.Padding = [0 5.5 0 5.5];
    presetHeaderCell.BackgroundColor = [1 1 1];
    ui.presetDropdown = uidropdown(presetHeaderCell, ...
        'Items', {'Custom', 'Square', 'Triangular', 'Honeycomb', 'Jx', 'Lieb', 'Ring'}, ...
        'Value', 'Custom', ...
        'BackgroundColor', [1 1 1], ...
        'Tag', 'PresetDropdown', ...
        'ValueChangedFcn', @presetChanged);
    ui.presetDropdown.Layout.Row = 1;
    ui.presetDropdown.Layout.Column = 1;

    sublatticePropertiesGrid = uigridlayout(sublatticePropertiesBody, [5 1]);
    sublatticePropertiesGrid.Layout.Row = 1;
    sublatticePropertiesGrid.Layout.Column = 1;
    sublatticePropertiesGrid.ColumnWidth = {'1x'};
    sublatticePropertiesGrid.RowHeight = {30, 26, 26, 26, 26};
    sublatticePropertiesGrid.Padding = [7 2 7 3];
    sublatticePropertiesGrid.RowSpacing = 3;
    sublatticePropertiesGrid.ColumnSpacing = 0;

    % The displayed equation is driven by the same profile state used by
    % automatic coupling evaluation and manual distance inversion.
    ui.formulaPanel = uipanel(sublatticePropertiesGrid, ...
        'BorderType', 'line', ...
        'BorderWidth', 1, ...
        'BorderColor', [0.86 0.86 0.86], ...
        'BackgroundColor', [0.97 0.975 0.98], ...
        'Tag', 'GeometricLawPanel');
    ui.formulaPanel.Layout.Row = 1;
    formulaGrid = uigridlayout(ui.formulaPanel, [1 1]);
    formulaGrid.Padding = [0 0 0 0];
    formulaGrid.RowHeight = {'1x'};
    formulaGrid.ColumnWidth = {'1x'};
    ui.formulaLabel = uilabel(formulaGrid, ...
        'Text', ['$\kappa_{i,j}=\exp(-\gamma r_{i,j})' ...
                 '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'], ...
        'Interpreter', 'latex', ...
        'FontSize', regularFontSize + 1, ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', [0.97 0.975 0.98], ...
        'Tag', 'GeometricLawLabel', ...
        'FontColor', [0.22 0.24 0.28]);

    profileRowGrid = uigridlayout(sublatticePropertiesGrid, [1 2]);
    profileRowGrid.Layout.Row = 2;
    profileRowGrid.ColumnWidth = {'1x', '1x'};
    profileRowGrid.RowHeight = {'1x'};
    profileRowGrid.Padding = [0 0 0 0];
    profileRowGrid.ColumnSpacing = 5;
    profileLabel = uilabel(profileRowGrid, 'Text', 'Coupling profile');
    profileLabel.Layout.Column = 1;
    ui.profileDropdown = uidropdown(profileRowGrid, ...
        'Items', {'Power law', 'Exponential', 'Fiber asymptotic'}, ...
        'Value', 'Exponential', ...
        'Tag', 'CouplingProfileDropdown', ...
        'ValueChangedFcn', @profileChanged);
    ui.profileDropdown.Layout.Column = 2;

    cutoffRowGrid = uigridlayout(sublatticePropertiesGrid, [1 2]);
    cutoffRowGrid.Layout.Row = 3;
    cutoffRowGrid.ColumnWidth = {'1x', '1x'};
    cutoffRowGrid.RowHeight = {'1x'};
    cutoffRowGrid.Padding = [0 0 0 0];
    cutoffRowGrid.ColumnSpacing = 5;

    gammaRowGrid = uigridlayout(sublatticePropertiesGrid, [1 2]);
    gammaRowGrid.Layout.Row = 4;
    gammaRowGrid.ColumnWidth = {'1x', '1x'};
    gammaRowGrid.RowHeight = {'1x'};
    gammaRowGrid.Padding = [0 0 0 0];
    gammaRowGrid.ColumnSpacing = 5;

    gammaLabel = uilabel(gammaRowGrid, ...
        'Text', 'Decay parameter, γ');
    gammaLabel.Layout.Row = 1;
    gammaLabel.Layout.Column = 1;
    ui.gammaField = uieditfield(gammaRowGrid, 'numeric', ...
        'Value', state.gamma, ...
        'Limits', [eps Inf], ...
        'Tag', 'GammaField', ...
        'Tooltip', 'Controls the distance decay of the selected coupling profile', ...
        'ValueChangedFcn', @gammaChanged);
    ui.gammaField.Layout.Row = 1;
    ui.gammaField.Layout.Column = 2;

    cutoffTooltip = ['Maximum absolute separation for automatic geometric ' ...
        'couplings; Inf disables truncation. Values below the shortest ' ...
        'edge are raised automatically.'];
    cutoffLabel = uilabel(cutoffRowGrid, ...
        'Text', 'Range cutoff', ...
        'Tooltip', cutoffTooltip);
    cutoffLabel.Layout.Row = 1;
    cutoffLabel.Layout.Column = 1;
    ui.cutoffField = uieditfield(cutoffRowGrid, 'numeric', ...
        'Value', state.cutoffDistance, ...
        'Tag', 'CutoffField', ...
        'Tooltip', cutoffTooltip, ...
        'ValueChangedFcn', @cutoffChanged);
    ui.cutoffField.Layout.Row = 1;
    ui.cutoffField.Layout.Column = 2;

    disorderRowGrid = uigridlayout(sublatticePropertiesGrid, [1 2]);
    disorderRowGrid.Layout.Row = 5;
    disorderRowGrid.ColumnWidth = {'1x', '1x'};
    disorderRowGrid.RowHeight = {'1x'};
    disorderRowGrid.Padding = [0 0 0 0];
    disorderRowGrid.ColumnSpacing = 5;
    disorderTypeGrid = uigridlayout(disorderRowGrid, [1 1]);
    disorderTypeGrid.Layout.Column = 1;
    disorderTypeGrid.ColumnWidth = {'1x'};
    disorderTypeGrid.RowHeight = {'1x'};
    disorderTypeGrid.Padding = [0 0 0 0];
    disorderTypeGrid.ColumnSpacing = 0;
    ui.disorderDropdown = uidropdown(disorderTypeGrid, ...
        'Items', cellstr(disorderTypes()), ...
        'Value', char(state.disorder.type), ...
        'Tag', 'DisorderTypeDropdown', ...
        'Tooltip', 'Select the disorder model', ...
        'ValueChangedFcn', @disorderTypeChanged);
    ui.disorderDropdown.Layout.Row = 1;
    ui.disorderDropdown.Layout.Column = 1;
    ui.disorderStrengthField = uieditfield(disorderRowGrid, 'numeric', ...
        'Value', state.disorder.strength, ...
        'Limits', [0 Inf], ...
        'Editable', 'off', ...
        'Enable', 'off', ...
        'BackgroundColor', readOnlyBackgroundColor, ...
        'FontColor', readOnlyFontColor, ...
        'Tag', 'DisorderStrengthField', ...
        'Tooltip', 'Disorder strength', ...
        'ValueChangedFcn', @disorderStrengthChanged);
    ui.disorderStrengthField.Layout.Row = 1;
    ui.disorderStrengthField.Layout.Column = 2;

    %% Element and preset-specific property controls

    % Custom uses this panel for node/edge editing. A named preset shows its
    % family parameters only while no element is selected; selection exposes
    % the same node/edge fields without converting the underlying preset.
    [selectionPanel,selectionBody,ui.selectionTitleLabel] = ...
        createSectionPanel(leftGrid, 'Element properties', ...
        'ConditionalPropertiesPanel');
    selectionPanel.Layout.Row = 2;
    selectionRootGrid = uigridlayout(selectionBody, [1 1]);
    selectionRootGrid.RowHeight = {'1x'};
    selectionRootGrid.ColumnWidth = {'1x'};
    selectionRootGrid.Padding = [0 0 0 0];
    selectionGrid = uigridlayout(selectionRootGrid, [4 4]);
    selectionGrid.Layout.Row = 1;
    selectionGrid.Layout.Column = 1;
    selectionGrid.ColumnWidth = {'1x', 100, 100, 100};
    selectionGrid.RowHeight = {'1x', 0, 0, 0};
    selectionGrid.Padding = [7 2 7 2];
    selectionGrid.RowSpacing = 3;
    selectionGrid.ColumnSpacing = 5;

    presetPropertiesGrid = uigridlayout(selectionRootGrid, [2 1]);
    presetPropertiesGrid.Layout.Row = 1;
    presetPropertiesGrid.Layout.Column = 1;
    presetPropertiesGrid.RowHeight = {26, 26};
    presetPropertiesGrid.ColumnWidth = {'1x'};
    presetPropertiesGrid.Padding = [7 25 7 25];
    presetPropertiesGrid.RowSpacing = 15;
    presetPropertiesGrid.ColumnSpacing = 0;
    presetPropertiesGrid.Visible = 'off';
    presetPropertiesGrid.Tag = 'PresetPropertiesGrid';
    presetParameterRowGrids = gobjects(1, 2);
    for parameterRowIndex = 1:2
        parameterRowGrid = uigridlayout(presetPropertiesGrid, [1 7]);
        parameterRowGrid.Layout.Row = parameterRowIndex;
        parameterRowGrid.Layout.Column = 1;
        parameterRowGrid.RowHeight = {'1x'};
        parameterRowGrid.ColumnWidth = ...
            {'1x', presetParameterFieldWidth, '1x', 0, 0, 0, 0};
        parameterRowGrid.Padding = [0 0 0 0];
        parameterRowGrid.ColumnSpacing = 0;
        presetParameterRowGrids(parameterRowIndex) = parameterRowGrid;
    end
    presetParameterMiddleGrid = uigridlayout(presetPropertiesGrid, [3 7]);
    presetParameterMiddleGrid.Layout.Row = [1 2];
    presetParameterMiddleGrid.Layout.Column = 1;
    presetParameterMiddleGrid.RowHeight = {'1x', 26, '1x'};
    presetParameterMiddleGrid.ColumnWidth = ...
        {'4x', presetParameterFieldWidth, '1x', ...
        presetParameterFieldWidth, '1x', presetParameterFieldWidth, '4x'};
    presetParameterMiddleGrid.Padding = [0 0 0 0];
    presetParameterMiddleGrid.RowSpacing = 0;
    presetParameterMiddleGrid.ColumnSpacing = 0;
    presetParameterMiddleGrid.Visible = 'off';
    ui.presetParameterCells = gobjects(1, 6);
    ui.presetParameterLabels = gobjects(1, 6);
    ui.presetParameterFields = gobjects(1, 6);
    for parameterControlIndex = 1:6
        parameterCell = uigridlayout(presetParameterRowGrids(1), [1 2]);
        parameterCell.Layout.Row = 1;
        parameterCell.Layout.Column = [1 2];
        parameterCell.RowHeight = {'1x'};
        parameterCell.ColumnWidth = {'1x', presetParameterFieldWidth};
        parameterCell.Padding = [0 0 0 0];
        parameterCell.ColumnSpacing = 4;
        parameterCell.Visible = 'off';
        ui.presetParameterCells(parameterControlIndex) = parameterCell;
        parameterLabel = uilabel(parameterCell, ...
            'Text', '', ...
            'Interpreter', 'tex', ...
            'HorizontalAlignment', 'right');
        parameterLabel.Layout.Row = 1;
        parameterLabel.Layout.Column = 1;
        ui.presetParameterLabels(parameterControlIndex) = parameterLabel;
        parameterField = uieditfield(parameterCell, 'numeric', ...
            'Value', 1, ...
            'Tag', sprintf('PresetParameter%dField', parameterControlIndex), ...
            'ValueChangedFcn', @(src, event) presetParameterChanged( ...
                src, event, parameterControlIndex));
        parameterField.Layout.Row = 1;
        parameterField.Layout.Column = 2;
        ui.presetParameterFields(parameterControlIndex) = parameterField;
    end

    ui.selectionHint = uilabel(selectionGrid, ...
        'Text', 'Select a site or edge in the editing plane', ...
        'WordWrap', 'on', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'center', ...
        'FontColor', [0.42 0.45 0.50]);
    ui.selectionHint.Layout.Row = [1 4];
    ui.selectionHint.Layout.Column = [1 4];

    ui.detuningLabel = uilabel(selectionGrid, 'Text', 'Detuning', 'Visible', 'off');
    ui.detuningLabel.Layout.Row = 1;
    ui.detuningLabel.Layout.Column = [1 3];
    ui.detuningField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Visible', 'off', ...
        'Tag', 'SiteDetuningField', ...
        'ValueChangedFcn', @siteDetuningChanged);
    ui.detuningField.Layout.Row = 1;
    ui.detuningField.Layout.Column = 4;

    ui.xLabel = uilabel(selectionGrid, 'Text', 'X-coordinate', 'Visible', 'off');
    ui.xLabel.Layout.Row = 2;
    ui.xLabel.Layout.Column = [1 3];
    ui.xField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Limits', [-5 5], ...
        'Visible', 'off', ...
        'Tag', 'SiteXField', ...
        'ValueChangedFcn', @siteXChanged);
    ui.xField.Layout.Row = 2;
    ui.xField.Layout.Column = 4;

    ui.yLabel = uilabel(selectionGrid, 'Text', 'Y-coordinate', 'Visible', 'off');
    ui.yLabel.Layout.Row = 3;
    ui.yLabel.Layout.Column = [1 3];
    ui.yField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Limits', [-5 5], ...
        'Visible', 'off', ...
        'Tag', 'SiteYField', ...
        'ValueChangedFcn', @siteYChanged);
    ui.yField.Layout.Row = 3;
    ui.yField.Layout.Column = 4;

    ui.directionHeader = uilabel(selectionGrid, ...
        'Text', 'Direction', ...
        'HorizontalAlignment', 'left', ...
        'Visible', 'off');
    ui.directionHeader.Layout.Row = 1;
    ui.directionHeader.Layout.Column = 1;
    ui.distanceHeader = uilabel(selectionGrid, ...
        'Text', 'Distance', ...
        'HorizontalAlignment', 'center', ...
        'Visible', 'off', ...
        'Tooltip', 'Directional effective distance in sublattice-coordinate units');
    ui.distanceHeader.Layout.Row = 1;
    ui.distanceHeader.Layout.Column = 2;
    ui.magnitudeHeader = uilabel(selectionGrid, ...
        'Text', 'Coupling', ...
        'HorizontalAlignment', 'center', ...
        'Visible', 'off');
    ui.magnitudeHeader.Layout.Row = 1;
    ui.magnitudeHeader.Layout.Column = 3;
    ui.phaseHeader = uilabel(selectionGrid, ...
        'Text', 'Phase (°)', ...
        'HorizontalAlignment', 'center', ...
        'Visible', 'off');
    ui.phaseHeader.Layout.Row = 1;
    ui.phaseHeader.Layout.Column = 4;

    ui.forwardLabel = uilabel(selectionGrid, 'Text', '1 → 2', 'Visible', 'off');
    ui.forwardLabel.Layout.Row = 2;
    ui.forwardLabel.Layout.Column = 1;
    ui.forwardDistanceField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Visible', 'off', ...
        'Tag', 'ForwardDistanceField', ...
        'Tooltip', 'Directional effective distance in sublattice-coordinate units', ...
        'ValueChangedFcn', @(src, event) directionalDistanceChanged(src, event, true));
    ui.forwardDistanceField.Layout.Row = 2;
    ui.forwardDistanceField.Layout.Column = 2;
    ui.forwardMagnitudeField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Limits', [0 Inf], ...
        'Visible', 'off', ...
        'Tag', 'ForwardMagnitudeField', ...
        'ValueChangedFcn', @(src, event) directionalPolarChanged(src, event, true, true));
    ui.forwardMagnitudeField.Layout.Row = 2;
    ui.forwardMagnitudeField.Layout.Column = 3;
    ui.forwardPhaseField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Visible', 'off', ...
        'Tag', 'ForwardPhaseField', ...
        'ValueChangedFcn', @(src, event) directionalPolarChanged(src, event, true, false));
    ui.forwardPhaseField.Layout.Row = 2;
    ui.forwardPhaseField.Layout.Column = 4;

    ui.backwardLabel = uilabel(selectionGrid, 'Text', '1 ← 2', 'Visible', 'off');
    ui.backwardLabel.Layout.Row = 3;
    ui.backwardLabel.Layout.Column = 1;
    ui.backwardDistanceField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Visible', 'off', ...
        'Tag', 'BackwardDistanceField', ...
        'Tooltip', 'Directional effective distance in sublattice-coordinate units', ...
        'ValueChangedFcn', @(src, event) directionalDistanceChanged(src, event, false));
    ui.backwardDistanceField.Layout.Row = 3;
    ui.backwardDistanceField.Layout.Column = 2;
    ui.backwardMagnitudeField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Limits', [0 Inf], ...
        'Visible', 'off', ...
        'Tag', 'BackwardMagnitudeField', ...
        'ValueChangedFcn', @(src, event) directionalPolarChanged(src, event, false, true));
    ui.backwardMagnitudeField.Layout.Row = 3;
    ui.backwardMagnitudeField.Layout.Column = 3;
    ui.backwardPhaseField = uieditfield(selectionGrid, 'numeric', ...
        'Value', 0, ...
        'Visible', 'off', ...
        'Tag', 'BackwardPhaseField', ...
        'ValueChangedFcn', @(src, event) directionalPolarChanged(src, event, false, false));
    ui.backwardPhaseField.Layout.Row = 3;
    ui.backwardPhaseField.Layout.Column = 4;

    ui.restoreButton = uibutton(selectionGrid, 'push', ...
        'Text', 'Restore geometric couplings', ...
        'Visible', 'off', ...
        'Tag', 'RestoreGeometricButton', ...
        'Tooltip', ['Restore the selected edge, or every coupling attached ' ...
                    'to the selected node, to automatic geometric values'], ...
        'ButtonPushedFcn', @restoreGeometricCouplings);
    ui.restoreButton.Layout.Row = 4;
    ui.restoreButton.Layout.Column = [1 4];

    %% Instructions and actions

    [instructionPanel,instructionBody,instructionTitleLabel] = ...
        createSectionPanel(leftGrid, 'Controls', 'SublatticeControlsPanel');
    instructionPanel.Layout.Row = 3;
    instructionGrid = uigridlayout(instructionBody, [4 1]);
    instructionGrid.RowHeight = {22, 108, 22, '1x'};
    instructionGrid.Padding = [8 6 8 6];
    instructionGrid.RowSpacing = 4;
    editingPlaneHeading = uilabel(instructionGrid, ...
        'Text', 'Editing plane', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
    editingPlaneHeading.Layout.Row = 1;
    editingPlaneDescription = uihtml(instructionGrid, ...
        'HTMLSource', justifiedControlsHtml( ...
            ['‘Double-click’ blank space to add a site. ‘Left-click’ to ' ...
             'select, ‘Shift-click’ or ‘Box-select’ for multiple sites, then ' ...
             'drag to move. ‘Right-click’, ‘Delete’, or ‘Backspace’ remove ' ...
             'sites, and ‘Ctrl+Z’ undoes the previous change. In the editing ' ...
             'plane, ‘Scroll’ to zoom, use the ‘Hand icon’ to pan, and the ' ...
             '‘House icon’ to recenter.']));
    editingPlaneDescription.Layout.Row = 2;
    propertiesPanelHeading = uilabel(instructionGrid, ...
        'Text', 'Properties panel', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom');
    propertiesPanelHeading.Layout.Row = 3;
    propertiesPanelDescription = uihtml(instructionGrid, ...
        'HTMLSource', justifiedControlsHtml( ...
            ['Selecting a preset replaces the Custom geometry. Editing a ' ...
             'preset site or edge converts it to Custom. Coupling profile, ' ...
             'range cutoff, and decay parameter modify couplings globally. ' ...
             'Disorder remains an overlay: detuning disorder is additive and ' ...
             'coupling disorder is multiplicative. Selecting a site or edge ' ...
             'exposes the corresponding editable fields.']));
    propertiesPanelDescription.Layout.Row = 4;

    buttonGrid = uigridlayout(leftGrid, [1 2]);
    buttonGrid.Layout.Row = 4;
    buttonGrid.ColumnWidth = {'1x', '1x'};
    buttonGrid.RowHeight = {'1x'};
    buttonGrid.Padding = [0 0 0 0];
    buttonGrid.ColumnSpacing = 8;
    ui.cancelButton = uibutton(buttonGrid, 'push', ...
        'Text', 'Cancel', ...
        'Tag', 'CancelSublatticeButton', ...
        'FontSize', buttonFontSize, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [1.00 0.86 0.86], ...
        'ButtonPushedFcn', @cancelConfiguration);
    ui.acceptButton = uibutton(buttonGrid, 'push', ...
        'Text', 'Build', ...
        'Tag', 'BuildSublatticeButton', ...
        'FontSize', buttonFontSize, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.84 0.95 0.84], ...
        'ButtonPushedFcn', @acceptConfiguration);

    %% Editing plane and view controls

    editorPanel = uipanel(mainGrid, ...
        'BorderType', 'none', ...
        'BackgroundColor', [1 1 1], ...
        'AutoResizeChildren', 'off');
    editorPanel.Layout.Row = 1;
    editorPanel.Layout.Column = 2;

    ui.axes = uiaxes(editorPanel, ...
        'Tag', 'SublatticeAxes', ...
        'Box', 'on', ...
        'Color', [1 1 1], ...
        'XColor', [0.38 0.41 0.46], ...
        'YColor', [0.38 0.41 0.46], ...
        'GridColor', [0.82 0.84 0.87], ...
        'MinorGridColor', [0.90 0.91 0.93], ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'XLim', [-5 5], ...
        'YLim', [-5 5], ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'FontSize', regularFontSize, ...
        'ButtonDownFcn', @axesClicked);
    xlabel(ui.axes, 'X', 'FontSize', axisLabelFontSize, 'Rotation', 0);
    ylabel(ui.axes, 'Y', 'FontSize', axisLabelFontSize, 'Rotation', 0);
    title(ui.axes, 'Sublattice editing plane', ...
        'FontSize', axisTitleFontSize, 'FontWeight', 'normal');
    ui.axes.PlotBoxAspectRatio = [1 1 1];
    ui.axes.PlotBoxAspectRatioMode = 'manual';
    disableDefaultInteractivity(ui.axes);
    ui.axes.Toolbar.Visible = 'off';
    editorIconDirectory = fullfile(matlabroot, 'toolbox', 'matlab', ...
        'uitools', 'uifigureappjs', 'release', 'toolbox', 'matlab', ...
        'rich_text_component', 'src', 'js', 'rtc', 'images', ...
        'figurefloatingpalette');
    panDefaultIcon = fullfile(editorIconDirectory, 'pan_default.svg');
    panActiveIcon = fullfile(editorIconDirectory, 'pan_active.svg');
    homeIcon = fullfile(editorIconDirectory, 'restoreview_action.svg');
    ui.panButton = uibutton(editorPanel, 'state', ...
        'Text', '✋', ...
        'Tag', 'PanViewButton', ...
        'FontSize', buttonFontSize, ...
        'Tooltip', 'Pan editing plane', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'ValueChangedFcn', @panModeChanged);
    if isfile(panDefaultIcon)
        ui.panButton.Icon = panDefaultIcon;
        ui.panButton.Text = '';
    end
    ui.recenterButton = uibutton(editorPanel, 'push', ...
        'Text', '⌂', ...
        'Tag', 'RecenterViewButton', ...
        'FontSize', buttonFontSize, ...
        'FontWeight', 'bold', ...
        'Tooltip', 'Recenter editing plane', ...
        'BackgroundColor', [0.97 0.97 0.97], ...
        'ButtonPushedFcn', @recenterEditingPlane);
    if isfile(homeIcon)
        ui.recenterButton.Icon = homeIcon;
        ui.recenterButton.Text = '';
    end

    %% Initialize and enter the blocking editor loop

    % Lay out and render once while hidden so the first visible frame uses
    % final typography, control sizes, and graph geometry.
    applyTypography();

    editorPanel.SizeChangedFcn = @layoutEditor;
    drawnow;
    if ~isvalid(ui.figure)
        return;
    end
    layoutEditor();
    if ~isempty(seededConfiguration)
        installSeededConfiguration(seededConfiguration);
    end
    refreshAll();
    state.initialized = true;
    if state.cancelRequested
        delete(ui.figure);
        return;
    end
    ui.figure.Visible = 'on';
    drawnow;
    if isvalid(ui.figure)
        layoutEditor();
        drawnow;
    end

    if ~isvalid(ui.figure)
        return;
    end

    try
        uiwait(ui.figure);
    catch exception
        if ~isvalid(ui.figure)
            return;
        end
        rethrow(exception);
    end

    if isvalid(ui.figure)
        delete(ui.figure);
    end
    clear editorFigureCleanup;

    %% Typography and responsive layout

    function [panel,bodyHost,titleLabel,headerGrid] = ...
            createSectionPanel(parent, titleText, tag)
        % Build a fixed-height white header, separator, and flexible body.
        panel = uipanel(parent, ...
            'BackgroundColor', [1 1 1], ...
            'Tag', tag);
        sectionGrid = uigridlayout(panel, [3 1]);
        sectionGrid.RowHeight = {sectionHeaderHeight - 1, 1, '1x'};
        sectionGrid.ColumnWidth = {'1x'};
        sectionGrid.Padding = [0 0 0 0];
        sectionGrid.RowSpacing = 0;

        headerGrid = uigridlayout(sectionGrid, [1 2]);
        headerGrid.Layout.Row = 1;
        headerGrid.Layout.Column = 1;
        headerGrid.ColumnWidth = {'1x', '1x'};
        headerGrid.RowHeight = {'1x'};
        headerGrid.Padding = [7 0 7 0];
        headerGrid.ColumnSpacing = 5;
        headerGrid.BackgroundColor = [1 1 1];
        headerGrid.Tag = [tag 'Header'];
        titleLabel = uilabel(headerGrid, ...
            'Text', titleText, ...
            'FontSize', panelTitleFontSize, ...
            'VerticalAlignment', 'center');
        titleLabel.Layout.Row = 1;
        titleLabel.Layout.Column = [1 2];

        separator = uipanel(sectionGrid, ...
            'BorderType', 'none', ...
            'BackgroundColor', sectionSeparatorColor, ...
            'Tag', [tag 'Separator']);
        separator.Layout.Row = 2;
        separator.Layout.Column = 1;

        bodyHost = uigridlayout(sectionGrid, [1 1]);
        bodyHost.Layout.Row = 3;
        bodyHost.Layout.Column = 1;
        bodyHost.RowHeight = {'1x'};
        bodyHost.ColumnWidth = {'1x'};
        bodyHost.Padding = [0 0 0 0];
    end

    function applyTypography()
        % Apply one body font first, then override headings and specialized text.
        fontHandles = findall(ui.figure, '-property', 'FontSize');
        for fontIndex = 1:numel(fontHandles)
            fontHandles(fontIndex).FontSize = regularFontSize;
        end
        sublatticePropertiesLabel.FontSize = panelTitleFontSize;
        ui.selectionTitleLabel.FontSize = panelTitleFontSize;
        instructionTitleLabel.FontSize = panelTitleFontSize;
        ui.cancelButton.FontSize = buttonFontSize;
        ui.acceptButton.FontSize = buttonFontSize;
        ui.restoreButton.FontSize = regularFontSize;
        ui.panButton.FontSize = buttonFontSize;
        ui.recenterButton.FontSize = buttonFontSize;
        ui.formulaLabel.FontSize = regularFontSize + 1;
        ui.axes.FontSize = regularFontSize;
        ui.axes.XLabel.FontSize = axisLabelFontSize;
        ui.axes.YLabel.FontSize = axisLabelFontSize;
        ui.axes.Title.FontSize = axisTitleFontSize;
        ui.axes.Title.FontWeight = 'normal';
    end

    function layoutEditor(~, ~)
        % Fit the square plotting region and align its compact view buttons.
        if ~isvalid(editorPanel) || ~isvalid(ui.axes)
            return;
        end
        if state.interaction.mode == "pan_pending" || state.interaction.mode == "pan"
            resetInteraction();
        end
        panelPosition = editorPanel.Position;
        width = max(panelPosition(3), 20);
        height = max(panelPosition(4), 24);
        ui.axes.Position = [10 10 max(width - 20, 1) max(height - 24, 1)];
        % Wait for UIAxes to expose the new InnerPosition before locating
        % controls against the square data frame. This keeps the buttons
        % attached to the frame throughout horizontal window resizing.
        drawnow nocallbacks;
        if ~isvalid(editorPanel) || ~isvalid(ui.axes)
            return;
        end
        buttonSize = 26;
        buttonGap = 2;
        buttonFrameGap = 5;
        plotBox = axesPlotBoxPosition();
        plotRight = sum(plotBox([1 3]));
        plotTop = sum(plotBox([2 4]));
        pairWidth = 2*buttonSize + buttonGap;
        pairLeft = max(10, plotRight - pairWidth);
        panX = pairLeft;
        homeX = pairLeft + buttonSize + buttonGap;
        buttonY = max(10, min(height - buttonSize, plotTop + buttonFrameGap));
        ui.panButton.Position = [panX, buttonY, buttonSize, buttonSize];
        ui.recenterButton.Position = [homeX, buttonY, buttonSize, buttonSize];
    end

    function position = axesPlotBoxPosition()
        % Return the data-bearing rectangle after equal-aspect letterboxing.
        position = ui.axes.InnerPosition;
        xSpan = diff(ui.axes.XLim);
        ySpan = diff(ui.axes.YLim);
        if xSpan <= 0 || ySpan <= 0
            return;
        end
        dataRatio = xSpan/ySpan;
        pixelRatio = position(3)/max(position(4), 1);
        if pixelRatio > dataRatio
            plotWidth = position(4)*dataRatio;
            position(1) = position(1) + (position(3) - plotWidth)/2;
            position(3) = plotWidth;
        elseif pixelRatio < dataRatio
            plotHeight = position(3)/dataRatio;
            position(2) = position(2) + (position(4) - plotHeight)/2;
            position(4) = plotHeight;
        end
    end

    %% View navigation

    function zoomEditingPlane(~, event)
        % Zoom about the pointer without changing lattice data or undo history.
        if ~isvalid(ui.axes) || state.interaction.mode ~= "none"
            return;
        end
        hitObject = hittest(ui.figure);
        if isempty(hitObject)
            return;
        elseif isequal(hitObject, ui.axes)
            hitAxes = ui.axes;
        else
            hitAxes = ancestor(hitObject, 'axes');
        end
        if isempty(hitAxes) || ~isequal(hitAxes, ui.axes)
            return;
        end
        anchor = double(ui.axes.CurrentPoint(1, 1:2));
        oldXLimits = ui.axes.XLim;
        oldYLimits = ui.axes.YLim;
        if any(~isfinite(anchor)) ...
                || anchor(1) < oldXLimits(1) || anchor(1) > oldXLimits(2) ...
                || anchor(2) < oldYLimits(1) || anchor(2) > oldYLimits(2)
            return;
        end
        zoomSteps = double(event.VerticalScrollCount);
        if ~isfinite(zoomSteps) || zoomSteps == 0
            return;
        end
        currentSpan = diff(oldXLimits);
        requestedScale = (16/15)^zoomSteps;
        targetSpan = min(200, max(0.5, currentSpan*requestedScale));
        scale = targetSpan/currentSpan;
        newXLimits = anchor(1) + (oldXLimits - anchor(1))*scale;
        newYLimits = anchor(2) + (oldYLimits - anchor(2))*scale;
        ui.axes.XLim = newXLimits;
        ui.axes.YLim = newYLimits;
        ui.axes.XLimMode = 'manual';
        ui.axes.YLimMode = 'manual';
        ui.xField.Limits = [min(ui.xField.Limits(1), newXLimits(1)), ...
            max(ui.xField.Limits(2), newXLimits(2))];
        ui.yField.Limits = [min(ui.yField.Limits(1), newYLimits(1)), ...
            max(ui.yField.Limits(2), newYLimits(2))];
        redrawGraph();
        drawnow limitrate nocallbacks;
    end

    function recenterEditingPlane(~, ~)
        % Restore the family-specific home view while retaining the lattice.
        applyHomeAxesLimits();
        redrawGraph();
        drawnow limitrate nocallbacks;
    end

    function panModeChanged(src, ~)
        % Pan mode changes gesture routing only; it is not configuration state.
        resetInteraction();
        clearPlacementClickCandidate();
        if src.Value
            if isfile(panActiveIcon)
                src.Icon = panActiveIcon;
            end
        else
            if isfile(panDefaultIcon)
                src.Icon = panDefaultIcon;
            end
        end
    end

    %% Keyboard shortcuts and undo history

    function handleUndoShortcut(sourceFigure, event)
        % Keep text editing native, and route lattice shortcuts at figure level.
        key = string(event.Key);
        modifiers = string(event.Modifier);
        isUndo = strcmpi(key, "z") ...
            && any(strcmpi(modifiers, "control") | strcmpi(modifiers, "command")) ...
            && ~any(strcmpi(modifiers, "shift"));
        if isUndo
            undoLastConfigurationChange();
            return;
        end
        isDelete = any(strcmpi(key, ["delete", "backspace"]));
        if isDelete && isempty(modifiers) && isequal(sourceFigure, ui.figure) ...
                && state.interaction.mode == "none" ...
                && ~keyboardFocusIsEditable(sourceFigure)
            deleteSites(state.selectedSiteUids);
        end
    end

    function tf = keyboardFocusIsEditable(figureHandle)
        % Protect active numeric/text entry from global undo and delete keys.
        focusedObject = figureHandle.CurrentObject;
        tf = ~isempty(focusedObject) && isvalid(focusedObject) ...
            && isprop(focusedObject, 'Editable') ...
            && strcmpi(string(focusedObject.Editable), "on");
    end

    function recordUndoState()
        % Capture immediately before a validated semantic mutation.
        pushUndoSnapshot(captureUndoState());
    end

    function pushUndoSnapshot(snapshot)
        % Retain a bounded history of semantic sublattice states.
        if isempty(snapshot)
            return;
        end
        undoStack{end + 1, 1} = snapshot;
        if numel(undoStack) > undoLimit
            undoStack(1) = [];
        end
    end

    function snapshot = captureUndoState()
        % View position and selection are intentionally absent.
        snapshot = struct( ...
            'nodes', state.nodes, ...
            'couplings', state.couplings, ...
            'manual', state.manual, ...
            'manualPhase', state.manualPhase, ...
            'syntheticDistanceRatio', state.syntheticDistanceRatio, ...
            'nextUid', state.nextUid, ...
            'couplingProfile', state.couplingProfile, ...
            'gamma', state.gamma, ...
            'cutoffDistance', state.cutoffDistance, ...
            'cutoffTracksPresetDefault', state.cutoffTracksPresetDefault, ...
            'disorder', state.disorder, ...
            'preset', state.preset, ...
            'presetSettings', state.presetSettings, ...
            'presetGlobals', state.presetGlobals, ...
            'customSnapshot', state.customSnapshot, ...
            'viewHomeHalfWidth', viewHomeHalfWidth);
    end

    function undoLastConfigurationChange()
        % Restore model state and view ownership, but never undo camera motion.
        if isempty(undoStack)
            return;
        end
        snapshot = undoStack{end};
        undoStack(end) = [];
        previousPreset = state.preset;
        previousUids = reshape([state.nodes.uid], 1, []);
        resetInteraction();
        clearPlacementClickCandidate();
        state.nodes = snapshot.nodes;
        state.couplings = snapshot.couplings;
        state.manual = snapshot.manual;
        state.manualPhase = snapshot.manualPhase;
        state.syntheticDistanceRatio = snapshot.syntheticDistanceRatio;
        state.nextUid = snapshot.nextUid;
        state.couplingProfile = snapshot.couplingProfile;
        state.gamma = snapshot.gamma;
        state.cutoffDistance = snapshot.cutoffDistance;
        state.cutoffTracksPresetDefault = snapshot.cutoffTracksPresetDefault;
        state.disorder = snapshot.disorder;
        state.preset = snapshot.preset;
        state.presetSettings = snapshot.presetSettings;
        state.presetGlobals = snapshot.presetGlobals;
        state.customSnapshot = snapshot.customSnapshot;
        viewHomeHalfWidth = snapshot.viewHomeHalfWidth;
        restoredUids = reshape([state.nodes.uid], 1, []);
        if previousPreset ~= state.preset ...
                || ~isequal(sort(previousUids), sort(restoredUids))
            clearSelection();
        else
            state.selectedSiteUids = state.selectedSiteUids( ...
                ismember(state.selectedSiteUids, restoredUids));
            if numel(state.selectedEdgeUids) ~= 2 ...
                    || ~all(ismember(state.selectedEdgeUids, restoredUids))
                state.selectedEdgeUids = [];
            end
        end
        removeUnavailableEdgeSelection();
        refreshAll();
    end

    %% Embedded sublattice configuration

    function configuration = sublatticeConfigurationFromState()
        % Capture base/manual state; automatic and effective values are rebuilt.
        count = numel(state.nodes);
        positions = zeros(count, 2);
        if count > 0
            positions(:, 1) = [state.nodes.x].';
            positions(:, 2) = [state.nodes.y].';
        end
        configuration = struct();
        configuration.format = 'sublattice_configuration';
        configuration.nextUid = state.nextUid;
        configuration.sites = struct( ...
            'uid', reshape([state.nodes.uid], [], 1), ...
            'position', positions, ...
            'detuning', reshape([state.nodes.detuning], [], 1));
        configuration.coupling = struct();
        configuration.coupling.profile = char(state.couplingProfile);
        configuration.coupling.gamma = state.gamma;
        configuration.coupling.cutoffDistance = state.cutoffDistance;
        configuration.coupling.manualMask = logical(state.manual);
        configuration.coupling.manualValue = complex(zeros(count));
        configuration.coupling.manualValue(state.manual) = state.couplings(state.manual);
        configuration.coupling.manualPhase = nan(count);
        configuration.coupling.manualPhase(state.manual) = state.manualPhase(state.manual);
        configuration.coupling.syntheticDistanceRatio = nan(count);
        configuration.coupling.syntheticDistanceRatio(state.manual) = ...
            state.syntheticDistanceRatio(state.manual);
        configuration.disorder = state.disorder;
        configuration.disorder.type = char(state.disorder.type);
        if state.preset == "Custom"
            presetParameters = struct();
        else
            presetParameters = state.presetSettings.(char(state.preset));
        end
        configuration.preset = struct( ...
            'name', char(state.preset), ...
            'parameters', presetParameters, ...
            'cutoffTracksDefault', state.cutoffTracksPresetDefault);
    end

    function configuration = validateSublatticeConfiguration(configuration)
        % Canonicalize only after every scalar, vector, and matrix invariant holds.
        requireConfiguration(isscalar(configuration) && isstruct(configuration), ...
            'The sublattice configuration must be a scalar structure.');
        requireFields(configuration, ...
            {'format', 'nextUid', 'sites', ...
             'coupling', 'disorder', 'preset'});
        requireConfiguration(strcmp(string(configuration.format), ...
            "sublattice_configuration"), ...
            'The record is not a sublattice configuration.');

        requireConfiguration(isscalar(configuration.sites) && isstruct(configuration.sites), ...
            'The saved site data are invalid.');
        requireFields(configuration.sites, {'uid', 'position', 'detuning'});
        uid = configuration.sites.uid;
        position = configuration.sites.position;
        detuning = configuration.sites.detuning;
        requireConfiguration(isnumeric(uid) && isreal(uid) && iscolumn(uid) ...
            && all(isfinite(uid)) && all(uid >= 1) && all(uid == fix(uid)) ...
            && all(uid <= flintmax) && numel(unique(uid)) == numel(uid), ...
            'Saved site identifiers are invalid.');
        count = numel(uid);
        requireConfiguration(isnumeric(position) && isreal(position) ...
            && isequal(size(position), [count 2]) && all(isfinite(position), 'all'), ...
            'Saved site coordinates are invalid.');
        requireConfiguration(isnumeric(detuning) && isreal(detuning) ...
            && iscolumn(detuning) && numel(detuning) == count ...
            && all(isfinite(detuning)), ...
            'Saved site detunings are invalid.');
        if count > 1
            for firstIndex = 1:count - 1
                separations = vecnorm(position(firstIndex + 1:end, :) ...
                    - position(firstIndex, :), 2, 2);
                requireConfiguration(all(separations >= 1e-8), ...
                    'Two saved sites occupy the same position.');
            end
        end
        nextUid = configuration.nextUid;
        minimumNextUid = 1;
        if ~isempty(uid)
            minimumNextUid = max(uid) + 1;
        end
        requireConfiguration(isnumeric(nextUid) && isreal(nextUid) ...
            && isscalar(nextUid) && isfinite(nextUid) && nextUid == fix(nextUid) ...
            && nextUid >= minimumNextUid && nextUid <= flintmax, ...
            'The saved next-site identifier is invalid.');

        requireConfiguration(isscalar(configuration.coupling) ...
            && isstruct(configuration.coupling), ...
            'The saved coupling data are invalid.');
        requireFields(configuration.coupling, {'profile', 'gamma', 'cutoffDistance', ...
            'manualMask', 'manualValue', 'manualPhase', 'syntheticDistanceRatio'});
        profile = string(configuration.coupling.profile);
        requireConfiguration(isscalar(profile) && any(profile == ...
            ["Power law", "Exponential", "Fiber asymptotic"]), ...
            'The saved coupling profile is invalid.');
        gamma = configuration.coupling.gamma;
        requireConfiguration(isnumeric(gamma) && isreal(gamma) && isscalar(gamma) ...
            && isfinite(gamma) && gamma > 0, ...
            'The saved decay parameter is invalid.');
        cutoffDistance = configuration.coupling.cutoffDistance;
        requireConfiguration(isnumeric(cutoffDistance) && isreal(cutoffDistance) ...
            && isscalar(cutoffDistance) && ~isnan(cutoffDistance) ...
            && cutoffDistance >= 0, ...
            'The saved range cutoff is invalid.');

        manualMask = configuration.coupling.manualMask;
        manualValue = configuration.coupling.manualValue;
        manualPhase = configuration.coupling.manualPhase;
        matrixSize = [count count];
        requireConfiguration((islogical(manualMask) || isnumeric(manualMask)) ...
            && isequal(size(manualMask), matrixSize) ...
            && all(manualMask(:) == 0 | manualMask(:) == 1), ...
            'The saved manual-coupling mask is invalid.');
        manualMask = logical(manualMask);
        requireConfiguration(~any(diag(manualMask)), ...
            'Saved self-couplings cannot be manually overridden.');
        requireConfiguration(isnumeric(manualValue) ...
            && isequal(size(manualValue), matrixSize) ...
            && all(isfinite(real(manualValue)), 'all') ...
            && all(isfinite(imag(manualValue)), 'all') ...
            && all(manualValue(~manualMask) == 0, 'all'), ...
            'The saved manual coupling values are invalid.');
        requireConfiguration(isnumeric(manualPhase) && isreal(manualPhase) ...
            && isequal(size(manualPhase), matrixSize) ...
            && all(isfinite(manualPhase(manualMask))) ...
            && all(isnan(manualPhase(~manualMask))), ...
            'The saved manual coupling phases are invalid.');
        syntheticDistanceRatio = configuration.coupling.syntheticDistanceRatio;
        requireConfiguration(isnumeric(syntheticDistanceRatio) ...
            && isreal(syntheticDistanceRatio) ...
            && isequal(size(syntheticDistanceRatio), matrixSize) ...
            && all((isfinite(syntheticDistanceRatio(manualMask)) ...
                & syntheticDistanceRatio(manualMask) >= 0) ...
                | (isinf(syntheticDistanceRatio(manualMask)) ...
                & syntheticDistanceRatio(manualMask) > 0)) ...
            && all(isnan(syntheticDistanceRatio(~manualMask))), ...
            'The saved synthetic coupling distance ratios are invalid.');
        disorder = validateDisorderConfiguration(configuration.disorder, uid);
        requireConfiguration(isscalar(configuration.preset) && isstruct(configuration.preset), ...
            'The saved preset data are invalid.');
        requireFields(configuration.preset, ...
            {'name', 'parameters', 'cutoffTracksDefault'});
        preset = string(configuration.preset.name);
        requireConfiguration(isscalar(preset), ...
            'The saved sublattice preset is invalid.');
        requireConfiguration(any(preset == ...
            ["Custom", "Square", "Triangular", "Honeycomb", "Jx", "Lieb", "Ring"]), ...
            'The saved sublattice preset is invalid.');
        presetParameters = validatePresetParameters( ...
            preset, configuration.preset.parameters);
        cutoffTracksDefault = configuration.preset.cutoffTracksDefault;
        requireConfiguration(islogical(cutoffTracksDefault) ...
            && isscalar(cutoffTracksDefault), ...
            'The saved preset-cutoff ownership is invalid.');
        requireConfiguration(preset ~= "Custom" || ~cutoffTracksDefault, ...
            'A Custom lattice cannot track a preset cutoff.');

        configuration.nextUid = double(nextUid);
        configuration.sites.uid = double(uid);
        configuration.sites.position = double(position);
        configuration.sites.detuning = double(detuning);
        configuration.coupling.profile = char(profile);
        configuration.coupling.gamma = double(gamma);
        configuration.coupling.cutoffDistance = double(cutoffDistance);
        configuration.coupling.manualMask = manualMask;
        configuration.coupling.manualValue = double(manualValue);
        configuration.coupling.manualPhase = double(manualPhase);
        configuration.coupling.syntheticDistanceRatio = double(syntheticDistanceRatio);
        configuration.disorder = disorder;
        configuration.preset.name = char(preset);
        configuration.preset.parameters = presetParameters;
        configuration.preset.cutoffTracksDefault = cutoffTracksDefault;
    end

    function disorder = validateDisorderConfiguration(disorder, activeUids)
        % Canonicalize one realization and verify UID/sample ownership.
        requireConfiguration(isscalar(disorder) && isstruct(disorder), ...
            'The saved disorder data are invalid.');
        requireFields(disorder, {'type', 'strength', 'siteUids', ...
            'siteSamples', 'edgeUidPairs', 'edgeSamples'});

        disorderType = string(disorder.type);
        requireConfiguration(isscalar(disorderType) ...
            && any(disorderType == disorderTypes()), ...
            'The saved disorder type is invalid.');
        strength = disorder.strength;
        requireConfiguration(isnumeric(strength) && isreal(strength) ...
            && isscalar(strength) && isfinite(strength) && strength >= 0, ...
            'The saved disorder strength is invalid.');
        requireConfiguration(disorderType ~= "No disorder" || strength == 0, ...
            'No disorder must have zero strength.');

        siteUids = disorder.siteUids;
        siteSamples = disorder.siteSamples;
        requireConfiguration(isnumeric(siteUids) && isreal(siteUids) ...
            && iscolumn(siteUids) && all(isfinite(siteUids)) ...
            && all(siteUids == fix(siteUids)) ...
            && numel(unique(siteUids)) == numel(siteUids) ...
            && all(ismember(siteUids, activeUids)), ...
            'The saved disorder site identifiers are invalid.');
        requireConfiguration(isnumeric(siteSamples) && isreal(siteSamples) ...
            && iscolumn(siteSamples) && numel(siteSamples) == numel(siteUids) ...
            && all(isfinite(siteSamples)) ...
            && all(siteSamples >= -0.5 & siteSamples <= 0.5), ...
            'The saved disorder site samples are invalid.');

        edgeUidPairs = disorder.edgeUidPairs;
        edgeSamples = disorder.edgeSamples;
        requireConfiguration(isnumeric(edgeUidPairs) && isreal(edgeUidPairs) ...
            && size(edgeUidPairs, 2) == 2 && all(isfinite(edgeUidPairs), 'all') ...
            && all(edgeUidPairs == fix(edgeUidPairs), 'all') ...
            && all(edgeUidPairs(:, 1) < edgeUidPairs(:, 2)) ...
            && all(ismember(edgeUidPairs, activeUids), 'all') ...
            && size(unique(edgeUidPairs, 'rows'), 1) == size(edgeUidPairs, 1), ...
            'The saved disorder edge identifiers are invalid.');
        requireConfiguration(isnumeric(edgeSamples) && isreal(edgeSamples) ...
            && iscolumn(edgeSamples) && numel(edgeSamples) == size(edgeUidPairs, 1) ...
            && all(isfinite(edgeSamples)) ...
            && all(edgeSamples >= -0.5 & edgeSamples <= 0.5), ...
            'The saved disorder edge samples are invalid.');

        switch disorderType
            case "Anderson disorder"
                requireConfiguration(isequal(sort(siteUids), sort(activeUids)) ...
                    && isempty(edgeUidPairs) && isempty(edgeSamples), ...
                    'The saved Anderson-disorder realization is incomplete.');
            case "Coupling disorder"
                expectedPairs = allDisorderUidPairs(activeUids);
                requireConfiguration(isempty(siteUids) && isempty(siteSamples) ...
                    && isequal(sortrows(edgeUidPairs), expectedPairs), ...
                    'The saved coupling-disorder realization is incomplete.');
            otherwise
                requireConfiguration(isempty(siteUids) && isempty(siteSamples) ...
                    && isempty(edgeUidPairs) && isempty(edgeSamples), ...
                    'The saved disorder realization is inconsistent with its type.');
        end

        disorder = struct( ...
            'type', disorderType, ...
            'strength', double(strength), ...
            'siteUids', double(siteUids), ...
            'siteSamples', double(siteSamples), ...
            'edgeUidPairs', double(edgeUidPairs), ...
            'edgeSamples', double(edgeSamples));
    end

    function parameters = validatePresetParameters(preset, parameters)
        % Require exactly the parameter set declared by the selected family.
        requireConfiguration(isscalar(parameters) && isstruct(parameters), ...
            'The saved preset parameters are invalid.');
        [keys, ~, ~, integerMask, lowerBounds] = ...
            presetParameterDefinition(preset);
        if preset == "Custom"
            requireConfiguration(isempty(fieldnames(parameters)), ...
                'Custom lattices cannot contain preset parameters.');
            parameters = struct();
            return;
        end
        requireConfiguration(numel(fieldnames(parameters)) == numel(keys), ...
            'The saved preset parameters do not match the selected preset.');
        canonicalParameters = struct();
        for parameterIndex = 1:numel(keys)
            key = keys{parameterIndex};
            requireConfiguration(isfield(parameters, key), ...
                sprintf('The saved preset is missing the "%s" parameter.', key));
            value = parameters.(key);
            requireConfiguration(isnumeric(value) && isreal(value) ...
                && isscalar(value) && isfinite(value) ...
                && value >= lowerBounds(parameterIndex), ...
                sprintf('The saved "%s" preset parameter is invalid.', key));
            if integerMask(parameterIndex)
                requireConfiguration(value == fix(value), ...
                    sprintf('The saved "%s" preset parameter must be an integer.', key));
            end
            canonicalParameters.(key) = double(value);
        end
        parameters = canonicalParameters;
    end

    function requireFields(value, names)
        % Exact field sets catch misspellings instead of silently ignoring them.
        requireConfiguration(numel(fieldnames(value)) == numel(names), ...
            'The sublattice configuration contains unexpected fields.');
        for nameIndex = 1:numel(names)
            requireConfiguration(isfield(value, names{nameIndex}), ...
                sprintf('The sublattice configuration is missing "%s".', ...
                names{nameIndex}));
        end
    end

    function requireConfiguration(condition, message)
        % Use one stable identifier for every embedded-schema failure.
        if ~condition
            error('lattice_config:InvalidSublatticeConfiguration', ...
                '%s', message);
        end
    end

    %% Preset and disorder defaults

    function settings = defaultPresetSettings()
        % Family parameters are independent so switching presets restores each set.
        settings = struct();
        settings.Square = struct('nx', 5, 'ny', 5, 'dMin', 1);
        settings.Triangular = struct('layers', 3, 'dMin', 1);
        settings.Honeycomb = struct( ...
            'layers', 2, 'dMin', 1, 'deltaA', 0, 'deltaB', 0);
        settings.Jx = struct('nx', 5, 'ny', 5, 'dMin', 1);
        settings.Lieb = struct( ...
            'nx', 3, 'ny', 3, 'dMin', 1, ...
            'deltaA', 0, 'deltaB', 0, 'deltaC', 0);
        settings.Ring = struct('n', 20, 'dMin', 1);
    end

    function values = disorderTypes()
        % This ordering is shared by dropdown presentation and validation.
        values = ["No disorder", "Anderson disorder", ...
            "Coupling disorder", "Aubry–André disorder"];
    end

    function disorder = defaultDisorderState()
        % Sample banks use stable UIDs rather than changeable display indices.
        disorder = struct( ...
            'type', "No disorder", ...
            'strength', 0, ...
            'siteUids', zeros(0, 1), ...
            'siteSamples', zeros(0, 1), ...
            'edgeUidPairs', zeros(0, 2), ...
            'edgeSamples', zeros(0, 1));
    end

    function globals = defaultPresetGlobals()
        % Absolute defaults include the intended shell but exclude the next one.
        defaultSettings = defaultPresetSettings();
        defaultValues = struct( ...
            'couplingProfile', "Exponential", ...
            'gamma', 1, ...
            'cutoffDistance', 1, ...
            'cutoffTracksDefault', true, ...
            'disorder', defaultDisorderState());
        globals = struct();
        presetNames = {'Square', 'Triangular', 'Honeycomb', 'Jx', 'Lieb', 'Ring'};
        for presetIndex = 1:numel(presetNames)
            presetName = presetNames{presetIndex};
            values = defaultValues;
            values.cutoffDistance = presetDefaultCutoffDistance( ...
                presetName, defaultSettings.(presetName));
            globals.(presetName) = values;
        end
    end

    function cutoffDistance = presetDefaultCutoffDistance(preset, settings)
        % Scale shell-safe preset ratios into the lattice's coordinate units.
        switch string(preset)
            case {"Square", "Jx", "Lieb"}
                cutoffRatio = 1.4;
            case {"Triangular", "Honeycomb"}
                cutoffRatio = 1.7;
            case "Ring"
                cutoffRatio = ringDefaultCutoffRatio(settings.n);
            otherwise
                error('lattice_config:UnknownPreset', ...
                    'Unknown lattice preset "%s".', preset);
        end
        cutoffDistance = cutoffRatio*settings.dMin;
    end

    function cutoffRatio = ringDefaultCutoffRatio(siteCount)
        % The first non-neighbor chord has ratio 2*cos(pi/N) when N exceeds 3.
        if siteCount <= 3
            cutoffRatio = 1;
            return;
        end
        firstUnwantedShell = 2*cos(pi/siteCount);
        cutoffRatio = max(1, floor(10*(firstUnwantedShell - 1e-10))/10);
    end

    function [keys, labels, tags, integerMask, lowerBounds] = ...
            presetParameterDefinition(preset)
        % Central metadata drives construction, display, callbacks, and validation.
        switch string(preset)
            case "Square"
                keys = {'nx', 'ny', 'dMin'};
                labels = {'N_x', 'N_y', 'd_{min}'};
                tags = {'SquareNxField', 'SquareNyField', 'SquareDMinField'};
                integerMask = [true true false];
                lowerBounds = [1 1 1e-6];
            case "Triangular"
                keys = {'layers', 'dMin'};
                labels = {'L', 'd_{min}'};
                tags = {'TriangularLayersField', 'TriangularDMinField'};
                integerMask = [true false];
                lowerBounds = [1 1e-6];
            case "Honeycomb"
                keys = {'layers', 'dMin', 'deltaA', 'deltaB'};
                labels = {'L', 'd_{min}', '\delta_A', '\delta_B'};
                tags = {'HoneycombLayersField', 'HoneycombDMinField', ...
                    'HoneycombDeltaAField', 'HoneycombDeltaBField'};
                integerMask = [true false false false];
                lowerBounds = [1 1e-6 -Inf -Inf];
            case "Jx"
                keys = {'nx', 'ny', 'dMin'};
                labels = {'N_x', 'N_y', 'd_{min}'};
                tags = {'JxNxField', 'JxNyField', 'JxDMinField'};
                integerMask = [true true false];
                lowerBounds = [1 1 1e-6];
            case "Lieb"
                keys = {'nx', 'ny', 'dMin', 'deltaA', 'deltaB', 'deltaC'};
                labels = {'N_x', 'N_y', 'd_{min}', '\delta_A', '\delta_B', '\delta_C'};
                tags = {'LiebNxField', 'LiebNyField', 'LiebDMinField', ...
                    'LiebDeltaAField', 'LiebDeltaBField', 'LiebDeltaCField'};
                integerMask = [true true false false false false];
                lowerBounds = [1 1 1e-6 -Inf -Inf -Inf];
            case "Ring"
                keys = {'n', 'dMin'};
                labels = {'N', 'd_{min}'};
                tags = {'RingNField', 'RingDMinField'};
                integerMask = [true false];
                lowerBounds = [3 1e-6];
            otherwise
                keys = {};
                labels = {};
                tags = {};
                integerMask = false(1, 0);
                lowerBounds = zeros(1, 0);
        end
    end

    function rows = presetParameterLayout(preset)
        % Assign each family's fixed-width fields to two available rows.
        switch string(preset)
            case {"Square", "Jx"}
                rows = [1 1 2];
            case "Triangular"
                rows = [1 2];
            case "Honeycomb"
                rows = [1 1 2 2];
            case "Lieb"
                rows = [1 1 1 2 2 2];
            case "Ring"
                rows = [1 2];
            otherwise
                rows = zeros(1, 0);
        end
    end

    %% Graph snapshots and preset lifecycle

    function snapshot = captureGraphSnapshot()
        % The Custom snapshot is a preset-switch cache, not an undo entry.
        snapshot = struct( ...
            'nodes', state.nodes, ...
            'couplings', state.couplings, ...
            'manual', state.manual, ...
            'manualPhase', state.manualPhase, ...
            'syntheticDistanceRatio', state.syntheticDistanceRatio, ...
            'nextUid', state.nextUid, ...
            'couplingProfile', state.couplingProfile, ...
            'gamma', state.gamma, ...
            'cutoffDistance', state.cutoffDistance, ...
            'cutoffTracksPresetDefault', state.cutoffTracksPresetDefault, ...
            'disorder', state.disorder);
    end

    function restoreGraphSnapshot(snapshot)
        % Reinstall a complete graph snapshot without changing persistence UI.
        state.nodes = snapshot.nodes;
        state.couplings = snapshot.couplings;
        state.manual = snapshot.manual;
        state.manualPhase = snapshot.manualPhase;
        state.syntheticDistanceRatio = snapshot.syntheticDistanceRatio;
        state.nextUid = snapshot.nextUid;
        state.couplingProfile = snapshot.couplingProfile;
        state.gamma = snapshot.gamma;
        state.cutoffDistance = snapshot.cutoffDistance;
        state.cutoffTracksPresetDefault = snapshot.cutoffTracksPresetDefault;
        state.disorder = snapshot.disorder;
        updateAxesLimits();
    end

    function resetToEmptyCustom()
        % Establish invariant empty matrices and the default Custom model.
        state.nodes = struct('uid', {}, 'x', {}, 'y', {}, 'detuning', {});
        state.couplings = complex(zeros(0));
        state.manual = false(0);
        state.manualPhase = nan(0);
        state.syntheticDistanceRatio = nan(0);
        state.nextUid = 1;
        state.couplingProfile = "Exponential";
        state.gamma = 1;
        state.cutoffDistance = Inf;
        state.cutoffTracksPresetDefault = false;
        state.disorder = defaultDisorderState();
        updateAxesLimits();
    end

    function cacheActivePresetGlobals()
        % Common settings are remembered independently for every preset family.
        if state.preset == "Custom"
            return;
        end
        state.presetGlobals.(char(state.preset)) = struct( ...
            'couplingProfile', state.couplingProfile, ...
            'gamma', state.gamma, ...
            'cutoffDistance', state.cutoffDistance, ...
            'cutoffTracksDefault', state.cutoffTracksPresetDefault, ...
            'disorder', state.disorder);
    end

    function activatePresetGlobals(preset)
        % Restore common coupling/disorder settings cached for this family.
        values = state.presetGlobals.(char(preset));
        state.couplingProfile = values.couplingProfile;
        state.gamma = values.gamma;
        state.cutoffDistance = values.cutoffDistance;
        state.cutoffTracksPresetDefault = values.cutoffTracksDefault;
        state.disorder = values.disorder;
    end

    function converted = convertPresetToCustomForEdit()
        % Convert in place and deliberately preserve the current camera limits.
        converted = state.preset ~= "Custom";
        if ~converted
            return;
        end
        cacheActivePresetGlobals();
        state.preset = "Custom";
        state.cutoffTracksPresetDefault = false;
        state.customSnapshot = [];
        ui.presetDropdown.Value = 'Custom';
    end

    function generateActivePreset(resetView)
        % Replace the graph from family parameters, then reconstruct auto edges.
        if nargin < 1
            resetView = true;
        end
        if state.preset == "Custom"
            return;
        end
        settings = state.presetSettings.(char(state.preset));
        [positions, detunings, engineeredEdges] = ...
            presetGeometry(state.preset, settings);
        count = size(positions, 1);
        generatedNodes = struct('uid', cell(1, count), 'x', cell(1, count), ...
            'y', cell(1, count), 'detuning', cell(1, count));
        for siteIndex = 1:count
            generatedNodes(siteIndex) = struct( ...
                'uid', siteIndex, ...
                'x', positions(siteIndex, 1), ...
                'y', positions(siteIndex, 2), ...
                'detuning', detunings(siteIndex));
        end
        state.nodes = generatedNodes;
        state.couplings = complex(zeros(count));
        state.manual = false(count);
        state.manualPhase = nan(count);
        state.syntheticDistanceRatio = nan(count);
        state.nextUid = count + 1;
        for edgeIndex = 1:size(engineeredEdges, 1)
            % Engineered Jx bonds are explicit symmetric overrides on a
            % uniformly spaced grid; ordinary pairs still follow geometry.
            firstIndex = engineeredEdges(edgeIndex, 1);
            secondIndex = engineeredEdges(edgeIndex, 2);
            magnitude = engineeredEdges(edgeIndex, 3);
            distanceRatio = inferSyntheticDistanceRatio(magnitude);
            state.couplings(firstIndex, secondIndex) = magnitude;
            state.couplings(secondIndex, firstIndex) = magnitude;
            state.manual(firstIndex, secondIndex) = true;
            state.manual(secondIndex, firstIndex) = true;
            state.manualPhase(firstIndex, secondIndex) = 0;
            state.manualPhase(secondIndex, firstIndex) = 0;
            state.syntheticDistanceRatio(firstIndex, secondIndex) = distanceRatio;
            state.syntheticDistanceRatio(secondIndex, firstIndex) = distanceRatio;
        end
        reorderSites();
        preserveCutoffDistance(state.cutoffDistance);
        recomputeAutomaticCouplings();
        reconcileDisorderRealization();
        cacheActivePresetGlobals();
        state.selectedSiteUids = [];
        state.selectedEdgeUids = [];
        clearPlacementClickCandidate();
        resetInteraction();
        if resetView
            updateAxesLimits();
        end
    end

    %% Disorder realization and effective lattice values

    function initializeDisorderRealization()
        % Entering a random family starts a fresh sample bank for that family.
        disorderType = state.disorder.type;
        disorderStrength = state.disorder.strength;
        state.disorder = defaultDisorderState();
        state.disorder.type = disorderType;
        state.disorder.strength = disorderStrength;
        reconcileDisorderRealization();
    end

    function reconcileDisorderRealization()
        % Preserve samples for surviving UIDs/pairs and draw only missing entries.
        uids = reshape([state.nodes.uid], [], 1);
        switch state.disorder.type
            case "Anderson disorder"
                samples = zeros(numel(uids), 1);
                [preserved, oldLocations] = ismember(uids, state.disorder.siteUids);
                samples(preserved) = state.disorder.siteSamples(oldLocations(preserved));
                samples(~preserved) = rand(nnz(~preserved), 1) - 0.5;
                state.disorder.siteUids = uids;
                state.disorder.siteSamples = samples;
                state.disorder.edgeUidPairs = zeros(0, 2);
                state.disorder.edgeSamples = zeros(0, 1);
            case "Coupling disorder"
                pairs = allDisorderUidPairs(uids);
                samples = zeros(size(pairs, 1), 1);
                [preserved, oldLocations] = ismember( ...
                    pairs, state.disorder.edgeUidPairs, 'rows');
                samples(preserved) = state.disorder.edgeSamples(oldLocations(preserved));
                samples(~preserved) = rand(nnz(~preserved), 1) - 0.5;
                state.disorder.siteUids = zeros(0, 1);
                state.disorder.siteSamples = zeros(0, 1);
                state.disorder.edgeUidPairs = pairs;
                state.disorder.edgeSamples = samples;
            otherwise
                state.disorder.siteUids = zeros(0, 1);
                state.disorder.siteSamples = zeros(0, 1);
                state.disorder.edgeUidPairs = zeros(0, 2);
                state.disorder.edgeSamples = zeros(0, 1);
        end
    end

    function pairs = allDisorderUidPairs(uids)
        % Canonical sorted UID pairs make samples independent of display order.
        uids = sort(reshape(uids, [], 1));
        if numel(uids) < 2
            pairs = zeros(0, 2);
        else
            pairs = nchoosek(uids, 2);
        end
    end

    function values = effectiveDetunings()
        % Apply onsite disorder lazily so node detunings stay editable base values.
        values = reshape([state.nodes.detuning], [], 1);
        switch state.disorder.type
            case "Anderson disorder"
                uids = reshape([state.nodes.uid], [], 1);
                [found, locations] = ismember(uids, state.disorder.siteUids);
                if ~all(found)
                    error('lattice_config:MissingDisorderSample', ...
                        'An Anderson-disorder site sample is missing.');
                end
                values = values + state.disorder.strength ...
                    * state.disorder.siteSamples(locations);
            case "Aubry–André disorder"
                beta = (sqrt(5) - 1)/2;
                indices = (1:numel(values)).';
                values = values + state.disorder.strength ...
                    * cos(2*pi*beta*indices);
        end
    end

    function values = effectiveCouplings()
        % One real factor is shared by both directions of each unordered pair.
        values = state.couplings;
        if state.disorder.type ~= "Coupling disorder"
            return;
        end
        pairs = state.disorder.edgeUidPairs;
        if isempty(pairs)
            return;
        end
        uids = reshape([state.nodes.uid], 1, []);
        [found, locations] = ismember(pairs, uids);
        if ~all(found, 'all')
            error('lattice_config:InvalidDisorderPair', ...
                'A coupling-disorder edge sample has invalid site identifiers.');
        end
        count = numel(uids);
        forwardIndices = sub2ind([count count], locations(:, 1), locations(:, 2));
        backwardIndices = sub2ind([count count], locations(:, 2), locations(:, 1));
        factors = 1 + state.disorder.strength*state.disorder.edgeSamples;
        values(forwardIndices) = factors.*values(forwardIndices);
        values(backwardIndices) = factors.*values(backwardIndices);
    end

    function ratio = displayedSyntheticDistanceRatio( ...
            sourceIndex, destinationIndex, effectiveValue)
        % Reflect coupling disorder in edge geometry without mutating base ratios.
        ratio = state.syntheticDistanceRatio(sourceIndex, destinationIndex);
        if state.disorder.type ~= "Coupling disorder"
            return;
        end

        baseMagnitude = abs(state.couplings(sourceIndex, destinationIndex));
        effectiveMagnitude = abs(effectiveValue);
        if isnan(ratio)
            ratio = inferSyntheticDistanceRatio(effectiveMagnitude);
            return;
        end
        if baseMagnitude == effectiveMagnitude || baseMagnitude == 0
            return;
        end

        baseProfileRatio = inferSyntheticDistanceRatio(baseMagnitude);
        effectiveProfileRatio = inferSyntheticDistanceRatio(effectiveMagnitude);
        if isinf(effectiveProfileRatio)
            ratio = Inf;
        elseif isfinite(baseProfileRatio) && isfinite(ratio)
            ratio = max(0, ratio + effectiveProfileRatio - baseProfileRatio);
        end
    end

    %% Preset geometry generators

    function [positions, detunings, engineeredEdges] = presetGeometry(preset, settings)
        % Return base coordinates/detunings and any explicitly engineered bonds.
        engineeredEdges = zeros(0, 3);
        switch string(preset)
            case "Square"
                positions = squarePositions(settings.nx, settings.ny, settings.dMin);
                detunings = zeros(size(positions, 1), 1);
            case "Triangular"
                % Axial coordinates clipped by a hexagonal layer radius.
                radius = settings.layers - 1;
                positions = zeros(1 + 3*settings.layers*(settings.layers - 1), 2);
                positionIndex = 0;
                for q = -radius:radius
                    for r = -radius:radius
                        if max([abs(q), abs(r), abs(q + r)]) <= radius
                            positionIndex = positionIndex + 1;
                            positions(positionIndex, :) = ...
                                [q + r/2, sqrt(3)*r/2]*settings.dMin;
                        end
                    end
                end
                positions = centeredPositions(positions(1:positionIndex, :));
                detunings = zeros(positionIndex, 1);
            case "Honeycomb"
                [positions, species] = honeycombPositions(settings.layers, settings.dMin);
                detunings = settings.deltaA*ones(size(positions, 1), 1);
                detunings(species == 2) = settings.deltaB;
            case "Jx"
                [positions, engineeredEdges] = ...
                    jxPositionsAndEdges(settings.nx, settings.ny, settings.dMin);
                detunings = zeros(size(positions, 1), 1);
            case "Lieb"
                [positions, detunings] = liebPositions(settings);
            case "Ring"
                % Choose the circumradius so adjacent chord length is d_min.
                siteCount = settings.n;
                radius = settings.dMin/(2*sin(pi/siteCount));
                angles = pi/2 - 2*pi*(0:siteCount - 1).'/siteCount;
                positions = radius*[cos(angles), sin(angles)];
                positions = centeredPositions(positions);
                detunings = zeros(siteCount, 1);
            otherwise
                error('lattice_config:UnknownPreset', ...
                    'Unknown lattice preset "%s".', preset);
        end
    end

    function positions = squarePositions(nx, ny, spacing)
        % Open Cartesian patch centered on its bounding box.
        [xIndex, yIndex] = ndgrid(0:nx - 1, 0:ny - 1);
        positions = [xIndex(:), -yIndex(:)]*spacing;
        positions = centeredPositions(positions);
    end

    function [positions, species] = honeycombPositions(layers, spacing)
        % Deduplicate exact integer vertex keys from fused hexagonal cells.
        radius = layers - 1;
        centers = zeros(1 + 3*layers*(layers - 1), 2);
        centerIndex = 0;
        for q = -radius:radius
            for r = -radius:radius
                if max([abs(q), abs(r), abs(q + r)]) <= radius
                    centerIndex = centerIndex + 1;
                    centers(centerIndex, :) = [2*q + r, 3*r];
                end
            end
        end
        centers = centers(1:centerIndex, :);
        vertexOffsets = [1 1; 0 2; -1 1; -1 -1; 0 -2; 1 -1];
        allKeys = zeros(6*centerIndex, 2);
        for centerNumber = 1:centerIndex
            rows = (centerNumber - 1)*6 + (1:6);
            allKeys(rows, :) = centers(centerNumber, :) + vertexOffsets;
        end
        [uniqueKeys, ~, keyMap] = unique(allKeys, 'rows');
        siteCount = size(uniqueKeys, 1);
        adjacency = false(siteCount);
        for centerNumber = 1:centerIndex
            rows = (centerNumber - 1)*6 + (1:6);
            vertexIndices = keyMap(rows);
            for vertexNumber = 1:6
                nextVertex = mod(vertexNumber, 6) + 1;
                adjacency(vertexIndices(vertexNumber), vertexIndices(nextVertex)) = true;
                adjacency(vertexIndices(nextVertex), vertexIndices(vertexNumber)) = true;
            end
        end
        positions = [sqrt(3)*uniqueKeys(:, 1)/2, uniqueKeys(:, 2)/2]*spacing;
        positions = centeredPositions(positions);
        species = zeros(siteCount, 1);
        % Breadth-first graph coloring assigns deterministic A/B sublattices.
        ordering = sortrows([(1:siteCount).', -positions(:, 2), positions(:, 1)], [2 3]);
        for orderedIndex = 1:siteCount
            startIndex = ordering(orderedIndex, 1);
            if species(startIndex) ~= 0
                continue;
            end
            species(startIndex) = 1;
            queue = startIndex;
            queueIndex = 1;
            while queueIndex <= numel(queue)
                currentIndex = queue(queueIndex);
                queueIndex = queueIndex + 1;
                neighbors = find(adjacency(currentIndex, :));
                for neighbor = neighbors
                    if species(neighbor) == 0
                        species(neighbor) = 3 - species(currentIndex);
                        queue(end + 1) = neighbor; %#ok<AGROW>
                    end
                end
            end
        end
    end

    function [positions, engineeredEdges] = jxPositionsAndEdges(nx, ny, spacing)
        % Normalize sqrt(n*(N-n)) axial weights by one global 2-D maximum.
        positions = squarePositions(nx, ny, spacing);
        indexGrid = reshape(1:nx*ny, [nx ny]);
        rawEdges = zeros(max(0, (nx - 1)*ny + nx*(ny - 1)), 3);
        edgeIndex = 0;
        for yIndex = 1:ny
            for xIndex = 1:nx - 1
                edgeIndex = edgeIndex + 1;
                rawEdges(edgeIndex, :) = [indexGrid(xIndex, yIndex), ...
                    indexGrid(xIndex + 1, yIndex), sqrt(xIndex*(nx - xIndex))];
            end
        end
        for yIndex = 1:ny - 1
            for xIndex = 1:nx
                edgeIndex = edgeIndex + 1;
                rawEdges(edgeIndex, :) = [indexGrid(xIndex, yIndex), ...
                    indexGrid(xIndex, yIndex + 1), sqrt(yIndex*(ny - yIndex))];
            end
        end
        rawEdges = rawEdges(1:edgeIndex, :);
        if isempty(rawEdges)
            engineeredEdges = zeros(0, 3);
        else
            rawEdges(:, 3) = rawEdges(:, 3)/max(rawEdges(:, 3));
            engineeredEdges = rawEdges;
        end
    end

    function [positions, detunings] = liebPositions(settings)
        % A sites use pitch 2*d_min; B and C occupy horizontal/vertical midpoints.
        nx = settings.nx;
        ny = settings.ny;
        positions = zeros(3*nx*ny - nx - ny, 2);
        detunings = zeros(size(positions, 1), 1);
        positionIndex = 0;
        for xIndex = 1:nx
            for yIndex = 1:ny
                positionIndex = positionIndex + 1;
                positions(positionIndex, :) = [2*(xIndex - 1), -2*(yIndex - 1)];
                detunings(positionIndex) = settings.deltaA;
            end
        end
        for xIndex = 1:nx - 1
            for yIndex = 1:ny
                positionIndex = positionIndex + 1;
                positions(positionIndex, :) = [2*xIndex - 1, -2*(yIndex - 1)];
                detunings(positionIndex) = settings.deltaB;
            end
        end
        for xIndex = 1:nx
            for yIndex = 1:ny - 1
                positionIndex = positionIndex + 1;
                positions(positionIndex, :) = [2*(xIndex - 1), -(2*yIndex - 1)];
                detunings(positionIndex) = settings.deltaC;
            end
        end
        positions = centeredPositions(positions*settings.dMin);
    end

    function positions = centeredPositions(positions)
        % Center a preset on the midpoint of its axis-aligned bounds.
        if isempty(positions)
            return;
        end
        boundsCenter = (min(positions, [], 1) + max(positions, [], 1))/2;
        positions = positions - boundsCenter;
    end

    %% View limits and property display formatting

    function updateAxesLimits()
        % Reset both the home width and visible limits for an explicit family view.
        viewHomeHalfWidth = presetViewHalfWidth();
        applyHomeAxesLimits();
    end

    function applyHomeAxesLimits()
        % Expand beyond the nominal family domain only when sites require it.
        baseHalfWidth = viewHomeHalfWidth;
        if isempty(state.nodes)
            halfWidth = baseHalfWidth;
        else
            positions = [[state.nodes.x].', [state.nodes.y].'];
            maximumExtent = max(abs(positions), [], 'all');
            if maximumExtent <= baseHalfWidth
                halfWidth = baseHalfWidth;
            else
                halfWidth = 1.12*maximumExtent + 0.5;
            end
        end
        ui.axes.XLim = [-halfWidth halfWidth];
        ui.axes.YLim = [-halfWidth halfWidth];
        ui.axes.XLimMode = 'manual';
        ui.axes.YLimMode = 'manual';
        ui.xField.Limits = [-halfWidth halfWidth];
        ui.yField.Limits = [-halfWidth halfWidth];
    end

    function halfWidth = presetViewHalfWidth()
        % Use a stable family-specific home domain independent of later edits.
        switch state.preset
            case {"Square", "Triangular", "Jx", "Lieb"}
                halfWidth = 3;
            case {"Honeycomb", "Ring"}
                halfWidth = 4;
            otherwise
                halfWidth = 5;
        end
    end

    %% State installation and property callbacks

    function installConfigurationModel(configuration)
        % Install validated base/manual data, then derive ordering and auto edges.
        count = numel(configuration.sites.uid);
        modelNodes = struct('uid', cell(1, count), 'x', cell(1, count), ...
            'y', cell(1, count), 'detuning', cell(1, count));
        for siteIndex = 1:count
            modelNodes(siteIndex) = struct( ...
                'uid', configuration.sites.uid(siteIndex), ...
                'x', configuration.sites.position(siteIndex, 1), ...
                'y', configuration.sites.position(siteIndex, 2), ...
                'detuning', configuration.sites.detuning(siteIndex));
        end
        state.nodes = modelNodes;
        state.couplings = configuration.coupling.manualValue;
        state.manual = configuration.coupling.manualMask;
        state.manualPhase = configuration.coupling.manualPhase;
        state.syntheticDistanceRatio = configuration.coupling.syntheticDistanceRatio;
        state.nextUid = configuration.nextUid;
        state.gamma = configuration.coupling.gamma;
        state.cutoffDistance = configuration.coupling.cutoffDistance;
        state.cutoffTracksPresetDefault = ...
            configuration.preset.cutoffTracksDefault;
        state.couplingProfile = string(configuration.coupling.profile);
        state.disorder = configuration.disorder;
        state.preset = string(configuration.preset.name);
        state.presetSettings = defaultPresetSettings();
        state.presetGlobals = defaultPresetGlobals();
        if state.preset ~= "Custom"
            state.presetSettings.(char(state.preset)) = ...
                configuration.preset.parameters;
        end
        state.customSnapshot = [];
        reorderSites();
        preserveCutoffDistance(state.cutoffDistance);
        cacheActivePresetGlobals();
        recomputeAutomaticCouplings();
    end

    function installSeededConfiguration(configuration)
        % Add GUI-only selection, interaction, and view cleanup after installation.
        installConfigurationModel(configuration);
        state.selectedSiteUids = [];
        state.selectedEdgeUids = [];
        clearPlacementClickCandidate();
        resetInteraction();
        removeUnavailableEdgeSelection();
        updateAxesLimits();
    end

    function presetChanged(src, ~)
        % Preserve the current Custom graph or preset globals before switching.
        newPreset = string(src.Value);
        if newPreset == state.preset
            return;
        end
        recordUndoState();
        oldPreset = state.preset;
        if oldPreset == "Custom"
            state.customSnapshot = captureGraphSnapshot();
        else
            cacheActivePresetGlobals();
        end
        if newPreset == "Custom"
            state.preset = "Custom";
            if isempty(state.customSnapshot)
                resetToEmptyCustom();
            else
                restoreGraphSnapshot(state.customSnapshot);
            end
            state.customSnapshot = [];
        else
            state.preset = newPreset;
            activatePresetGlobals(newPreset);
            generateActivePreset();
        end
        clearPlacementClickCandidate();
        resetInteraction();
        clearSelection();
        refreshAll();
    end

    function profileChanged(src, ~)
        % Presets regenerate engineered metadata; Custom retains manual entries.
        newProfile = string(src.Value);
        if newProfile == state.couplingProfile
            return;
        end
        recordUndoState();
        state.couplingProfile = newProfile;
        cacheActivePresetGlobals();
        if state.preset == "Custom"
            recomputeAutomaticCouplings();
            removeUnavailableEdgeSelection();
        else
            generateActivePreset(false);
        end
        refreshAll();
    end

    function gammaChanged(src, ~)
        % Validate the decay edit, then rebuild only profile-derived couplings.
        if ~isfinite(src.Value) || src.Value <= 0
            src.Value = state.gamma;
            return;
        end
        if src.Value == state.gamma
            return;
        end
        recordUndoState();
        state.gamma = src.Value;
        cacheActivePresetGlobals();
        if state.preset == "Custom"
            recomputeAutomaticCouplings();
            removeUnavailableEdgeSelection();
        else
            generateActivePreset(false);
        end
        refreshAll();
    end

    function cutoffChanged(src, ~)
        % Clamp the absolute cutoff to the current shortest-edge floor.
        value = src.Value;
        if isnan(value)
            value = state.cutoffDistance;
        else
            value = correctedCutoffDistance(max(0, value));
        end
        src.Value = value;
        if isequaln(value, state.cutoffDistance)
            return;
        end
        recordUndoState();
        state.cutoffDistance = value;
        state.cutoffTracksPresetDefault = false;
        cacheActivePresetGlobals();
        if state.preset == "Custom"
            recomputeAutomaticCouplings();
            removeUnavailableEdgeSelection();
        else
            generateActivePreset(false);
        end
        refreshAll();
    end

    function disorderTypeChanged(src, ~)
        % Changing disorder type creates a fresh realization on the base model.
        newType = string(src.Value);
        if newType == state.disorder.type
            return;
        end
        recordUndoState();
        state.disorder.type = newType;
        if newType == "No disorder"
            state.disorder.strength = 0;
        end
        initializeDisorderRealization();
        cacheActivePresetGlobals();
        refreshAll();
    end

    function disorderStrengthChanged(src, ~)
        % Change only the scale applied to the stored disorder realization.
        if state.disorder.type == "No disorder"
            src.Value = 0;
            return;
        end
        value = src.Value;
        if ~isfinite(value) || value < 0
            src.Value = state.disorder.strength;
            return;
        end
        if value == state.disorder.strength
            return;
        end
        recordUndoState();
        state.disorder.strength = value;
        cacheActivePresetGlobals();
        refreshAll();
    end

    function presetParameterChanged(src, ~, parameterIndex)
        % Validate one family parameter and regenerate the active preset.
        if state.preset == "Custom"
            return;
        end
        [keys, ~, ~, integerMask, lowerBounds] = ...
            presetParameterDefinition(state.preset);
        if parameterIndex > numel(keys)
            return;
        end
        settings = state.presetSettings.(char(state.preset));
        key = keys{parameterIndex};
        previousValue = settings.(key);
        value = src.Value;
        if ~isfinite(value)
            src.Value = previousValue;
            return;
        end
        if integerMask(parameterIndex)
            value = round(value);
        end
        if value < lowerBounds(parameterIndex)
            src.Value = previousValue;
            return;
        end
        src.Value = value;
        if isequaln(value, previousValue)
            return;
        end
        updatePresetCutoff = state.cutoffTracksPresetDefault ...
            && (strcmp(key, 'dMin') ...
            || (state.preset == "Ring" && strcmp(key, 'n')));
        % Untouched preset cutoffs track spacing and Ring shell geometry.
        recordUndoState();
        settings.(key) = value;
        state.presetSettings.(char(state.preset)) = settings;
        if updatePresetCutoff
            state.cutoffDistance = presetDefaultCutoffDistance( ...
                state.preset, settings);
            cacheActivePresetGlobals();
        end
        resetView = ~startsWith(string(key), "delta");
        generateActivePreset(resetView);
        refreshAll();
    end

    %% Editing-plane interaction state machine

    function axesClicked(~, event)
        % Blank left gestures become placement clicks or box selections.
        if ui.panButton.Value
            beginPan(event);
            return;
        end
        if mouseButton(event) ~= 1
            clearPlacementClickCandidate();
            return;
        end
        selectionType = string(ui.figure.SelectionType);
        position = clampedAxesPoint(clickPosition(event));
        if isPlacementDoubleClick(selectionType, position)
            resetInteraction();
            addSite(position);
            return;
        end
        if selectionType ~= "normal" && selectionType ~= "extend"
            return;
        end
        state.interaction.mode = "box_pending";
        state.interaction.startPoint = position;
        state.interaction.currentPoint = position;
        state.interaction.dragUids = [];
        state.interaction.originalPositions = zeros(0, 2);
        state.interaction.moved = false;
        state.interaction.additive = shiftIsHeld() || selectionType == "extend";
        if state.interaction.additive
            clearPlacementClickCandidate();
        end
        state.interaction.initialSelection = state.selectedSiteUids;
        state.interaction.clickedUid = [];
    end

    function nodeClicked(~, event, uid)
        % Delay mutation until movement exceeds the click/drag threshold.
        clearPlacementClickCandidate();
        button = mouseButton(event);
        if ui.panButton.Value && button == 1
            beginPan(event);
            return;
        end
        if button == 3
            resetInteraction();
            deleteSites(uid);
            return;
        elseif button ~= 1
            return;
        end
        selectionType = string(ui.figure.SelectionType);
        if selectionType == "open"
            return;
        elseif selectionType ~= "normal" && selectionType ~= "extend"
            return;
        end
        initialSelection = state.selectedSiteUids;
        isAdditive = shiftIsHeld() || selectionType == "extend";
        if isAdditive
            dragSelection = unique([initialSelection, uid], 'stable');
        elseif any(initialSelection == uid)
            dragSelection = initialSelection;
        else
            dragSelection = uid;
        end
        dragUndoSnapshot = captureUndoState();
        dragUndoRecorded = false;
        dragCutoffDistance = state.cutoffDistance;
        state.selectedSiteUids = dragSelection;
        state.selectedEdgeUids = [];
        state.interaction.mode = "node_pending";
        state.interaction.startPoint = clampedAxesPoint(clickPosition(event));
        state.interaction.currentPoint = state.interaction.startPoint;
        state.interaction.dragUids = dragSelection;
        state.interaction.originalPositions = positionsFromUids(dragSelection);
        state.interaction.moved = false;
        state.interaction.additive = isAdditive;
        state.interaction.initialSelection = initialSelection;
        state.interaction.clickedUid = uid;
        refreshSelectionEditor();
        redrawGraph();
    end

    function edgeClicked(~, event, uidA, uidB)
        % Edge selection is UID-based so reordering cannot invalidate it.
        clearPlacementClickCandidate();
        if ui.panButton.Value
            beginPan(event);
            return;
        end
        if mouseButton(event) ~= 1
            return;
        end
        resetInteraction();
        if selectedEdgeMatches(uidA, uidB)
            state.selectedEdgeUids = [];
        else
            state.selectedEdgeUids = [uidA uidB];
        end
        state.selectedSiteUids = [];
        refreshSelectionEditor();
        redrawGraph();
    end

    function beginPan(event)
        % Pixel-space deltas avoid feedback as axes limits move under the pointer.
        if mouseButton(event) ~= 1
            return;
        end
        resetInteraction();
        clearPlacementClickCandidate();
        state.interaction.mode = "pan_pending";
        state.interaction.startPoint = double(ui.figure.CurrentPoint(1, 1:2));
        state.interaction.currentPoint = state.interaction.startPoint;
        state.interaction.originalXLimits = ui.axes.XLim;
        state.interaction.originalYLimits = ui.axes.YLim;
        ui.figure.Pointer = 'hand';
    end

    function pointerMoved(~, ~)
        % Dispatch by latched mode; only the first real drag motion records undo.
        mode = state.interaction.mode;
        if mode == "none"
            return;
        end
        if mode == "pan_pending" || mode == "pan"
            currentPixel = double(ui.figure.CurrentPoint(1, 1:2));
            pixelMovement = currentPixel - state.interaction.startPoint;
            state.interaction.currentPoint = currentPixel;
            if ~state.interaction.moved && hypot(pixelMovement(1), pixelMovement(2)) < 1
                return;
            end
            state.interaction.mode = "pan";
            state.interaction.moved = true;
            plotBox = axesPlotBoxPosition();
            xShift = pixelMovement(1)*diff(state.interaction.originalXLimits) ...
                /max(plotBox(3), 1);
            yShift = pixelMovement(2)*diff(state.interaction.originalYLimits) ...
                /max(plotBox(4), 1);
            ui.axes.XLim = state.interaction.originalXLimits - xShift;
            ui.axes.YLim = state.interaction.originalYLimits - yShift;
            ui.axes.XLimMode = 'manual';
            ui.axes.YLimMode = 'manual';
            ui.xField.Limits = [min(ui.xField.Limits(1), ui.axes.XLim(1)), ...
                max(ui.xField.Limits(2), ui.axes.XLim(2))];
            ui.yField.Limits = [min(ui.yField.Limits(1), ui.axes.YLim(1)), ...
                max(ui.yField.Limits(2), ui.axes.YLim(2))];
            drawnow limitrate nocallbacks;
            return;
        end
        currentPoint = clampedAxesPoint(ui.axes.CurrentPoint(1, 1:2));
        state.interaction.currentPoint = currentPoint;
        movement = currentPoint - state.interaction.startPoint;
        if ~state.interaction.moved && ~movementExceedsThreshold(movement)
            return;
        end
        state.interaction.moved = true;
        state.lastBlankClickTimer = [];

        if mode == "box_pending" || mode == "box"
            state.interaction.mode = "box";
            redrawGraph();
            drawnow limitrate nocallbacks;
            return;
        end
        if mode ~= "node_pending" && mode ~= "node_drag"
            return;
        end
        state.interaction.mode = "node_drag";
        translatedPositions = clampedGroupPositions( ...
            state.interaction.originalPositions, movement);
        if positionsCollideWithUnselected(translatedPositions, state.interaction.dragUids)
            return;
        end
        if isequaln(translatedPositions, positionsFromUids(state.interaction.dragUids))
            return;
        end
        if ~dragUndoRecorded
            if isempty(dragUndoSnapshot)
                dragUndoSnapshot = captureUndoState();
            end
            pushUndoSnapshot(dragUndoSnapshot);
            dragUndoRecorded = true;
        end
        convertPresetToCustomForEdit();
        for dragIndex = 1:numel(state.interaction.dragUids)
            nodeIndex = siteIndexFromUid(state.interaction.dragUids(dragIndex));
            state.nodes(nodeIndex).x = translatedPositions(dragIndex, 1);
            state.nodes(nodeIndex).y = translatedPositions(dragIndex, 2);
        end
        preserveCutoffDistance(dragCutoffDistance);
        recomputeAutomaticCouplings();
        refreshPresetControls();
        refreshSelectionEditor();
        redrawGraph();
        drawnow limitrate nocallbacks;
    end

    function pointerReleased(~, ~)
        % Commit click/box/drag semantics, then return to the neutral mode.
        mode = state.interaction.mode;
        if mode == "none"
            return;
        end
        if mode == "pan_pending" || mode == "pan"
            resetInteraction();
            return;
        end
        interaction = state.interaction;
        if mode == "box_pending"
            if ~interaction.additive
                clearSelection();
            end
        elseif mode == "box"
            xBounds = sort([interaction.startPoint(1), interaction.currentPoint(1)]);
            yBounds = sort([interaction.startPoint(2), interaction.currentPoint(2)]);
            insideUids = [];
            for nodeIndex = 1:numel(state.nodes)
                if state.nodes(nodeIndex).x >= xBounds(1) && state.nodes(nodeIndex).x <= xBounds(2) ...
                        && state.nodes(nodeIndex).y >= yBounds(1) && state.nodes(nodeIndex).y <= yBounds(2)
                    insideUids(end + 1) = state.nodes(nodeIndex).uid; %#ok<AGROW>
                end
            end
            if interaction.additive
                state.selectedSiteUids = unique([interaction.initialSelection, insideUids], 'stable');
            else
                state.selectedSiteUids = insideUids;
            end
            state.selectedEdgeUids = [];
        elseif mode == "node_pending"
            if interaction.additive && any(interaction.initialSelection == interaction.clickedUid)
                state.selectedSiteUids = interaction.initialSelection( ...
                    interaction.initialSelection ~= interaction.clickedUid);
            elseif ~interaction.additive && isscalar(interaction.initialSelection) ...
                    && interaction.initialSelection == interaction.clickedUid
                state.selectedSiteUids = [];
            elseif ~interaction.additive && numel(interaction.initialSelection) > 1 ...
                    && any(interaction.initialSelection == interaction.clickedUid)
                state.selectedSiteUids = interaction.clickedUid;
            end
        elseif mode == "node_drag"
            reorderSites();
            recomputeAutomaticCouplings();
            removeUnavailableEdgeSelection();
        end
        resetInteraction();
        refreshAll();
    end

    function resetInteraction()
        % Clear transient gesture data without touching selection or view state.
        state.interaction.mode = "none";
        state.interaction.startPoint = [nan nan];
        state.interaction.currentPoint = [nan nan];
        state.interaction.dragUids = [];
        state.interaction.originalPositions = zeros(0, 2);
        state.interaction.originalXLimits = [nan nan];
        state.interaction.originalYLimits = [nan nan];
        state.interaction.moved = false;
        state.interaction.additive = false;
        state.interaction.initialSelection = [];
        state.interaction.clickedUid = [];
        dragUndoSnapshot = [];
        dragUndoRecorded = false;
        dragCutoffDistance = [];
        if isfield(ui, 'figure') && isvalid(ui.figure)
            ui.figure.Pointer = 'arrow';
        end
    end

    function clearPlacementClickCandidate()
        % Cancel the pending first click of a blank-space double-click.
        state.lastBlankClickTimer = [];
        state.lastBlankClickPosition = [nan nan];
    end

    function tf = isPlacementDoubleClick(selectionType, position)
        % Recognize a spatially local double-click without trusting event order.
        if isempty(state.lastBlankClickTimer)
            elapsedSeconds = Inf;
        else
            elapsedSeconds = toc(state.lastBlankClickTimer);
        end
        repeatedClick = elapsedSeconds <= 0.5 ...
            && ~movementExceedsThreshold(position - state.lastBlankClickPosition);
        tf = selectionType == "open" || ...
            (selectionType == "normal" && ~shiftIsHeld() && repeatedClick);
        if tf
            state.lastBlankClickTimer = [];
            state.lastBlankClickPosition = [nan nan];
        elseif selectionType == "normal" && ~shiftIsHeld()
            state.lastBlankClickTimer = tic;
            state.lastBlankClickPosition = position;
        end
    end

    function tf = shiftIsHeld()
        % Treat either Shift key as the additive-selection modifier.
        modifiers = ui.figure.CurrentModifier;
        tf = any(strcmpi(modifiers, 'shift'));
    end

    function position = clampedAxesPoint(position)
        % Keep edits inside the visible square domain.
        position = double(position(1:2));
        position(1) = min(ui.axes.XLim(2), max(ui.axes.XLim(1), position(1)));
        position(2) = min(ui.axes.YLim(2), max(ui.axes.YLim(1), position(2)));
    end

    function tf = movementExceedsThreshold(movement)
        % Convert the four-pixel gesture threshold into current data units.
        xPerPixel = diff(ui.axes.XLim)/max(ui.axes.Position(3), 1);
        yPerPixel = diff(ui.axes.YLim)/max(ui.axes.Position(4), 1);
        pixelMovement = hypot(movement(1)/xPerPixel, movement(2)/yPerPixel);
        tf = pixelMovement >= 4;
    end

    function positions = positionsFromUids(uids)
        % Resolve stable site identifiers to the current ordered coordinates.
        positions = zeros(numel(uids), 2);
        for uidIndex = 1:numel(uids)
            nodeIndex = siteIndexFromUid(uids(uidIndex));
            positions(uidIndex, :) = [state.nodes(nodeIndex).x, state.nodes(nodeIndex).y];
        end
    end

    function positions = clampedGroupPositions(originalPositions, movement)
        % Apply one translation while preserving group shape inside the axes.
        if isempty(originalPositions)
            positions = originalPositions;
            return;
        end
        movement(1) = min(ui.axes.XLim(2) - max(originalPositions(:, 1)), ...
            max(ui.axes.XLim(1) - min(originalPositions(:, 1)), movement(1)));
        movement(2) = min(ui.axes.YLim(2) - max(originalPositions(:, 2)), ...
            max(ui.axes.YLim(1) - min(originalPositions(:, 2)), movement(2)));
        positions = originalPositions + movement;
    end

    function tf = positionsCollideWithUnselected(candidatePositions, movingUids)
        % Reject drags that place a moved site on an unmoved site.
        tf = false;
        for nodeIndex = 1:numel(state.nodes)
            if any(movingUids == state.nodes(nodeIndex).uid)
                continue;
            end
            stationaryPosition = [state.nodes(nodeIndex).x, state.nodes(nodeIndex).y];
            if any(vecnorm(candidatePositions - stationaryPosition, 2, 2) < 1e-8)
                tf = true;
                return;
            end
        end
    end

    %% Site creation, deletion, and editing

    function addSite(position)
        % Extend every directional matrix and assign samples for new UID pairs.
        position = double(position(:).');
        if numel(position) ~= 2 || any(~isfinite(position))
            return;
        end
        position(1) = min(ui.axes.XLim(2), max(ui.axes.XLim(1), position(1)));
        position(2) = min(ui.axes.YLim(2), max(ui.axes.YLim(1), position(2)));
        if ~isempty(state.nodes)
            coordinates = [[state.nodes.x].', [state.nodes.y].'];
            if any(vecnorm(coordinates - position, 2, 2) < 1e-8)
                return;
            end
        end

        previousCutoffDistance = state.cutoffDistance;
        recordUndoState();
        convertPresetToCustomForEdit();
        uid = state.nextUid;
        state.nextUid = state.nextUid + 1;
        state.nodes(end + 1) = struct( ...
            'uid', uid, 'x', position(1), 'y', position(2), 'detuning', 0);
        oldCount = size(state.couplings, 1);
        newCount = oldCount + 1;
        state.couplings(newCount, newCount) = complex(0);
        state.manual(newCount, newCount) = false;
        state.manualPhase(newCount, newCount) = nan;
        state.manualPhase(newCount, :) = nan;
        state.manualPhase(:, newCount) = nan;
        state.syntheticDistanceRatio(newCount, newCount) = nan;
        state.syntheticDistanceRatio(newCount, :) = nan;
        state.syntheticDistanceRatio(:, newCount) = nan;
        reorderSites();
        reconcileDisorderRealization();
        preserveCutoffDistance(previousCutoffDistance);
        recomputeAutomaticCouplings();
        state.selectedSiteUids = uid;
        state.selectedEdgeUids = [];
        refreshAll();
    end

    function deleteSites(uids)
        % Delete a multi-selection as one mutation and one undo transaction.
        existingUids = reshape([state.nodes.uid], 1, []);
        uids = unique(reshape(uids, 1, []), 'stable');
        uids = uids(ismember(uids, existingUids));
        if isempty(uids)
            return;
        end
        previousCutoffDistance = state.cutoffDistance;
        recordUndoState();
        convertPresetToCustomForEdit();
        deleteMask = ismember([state.nodes.uid], uids);
        keepMask = ~deleteMask;
        state.nodes = state.nodes(keepMask);
        state.couplings = state.couplings(keepMask, keepMask);
        state.manual = state.manual(keepMask, keepMask);
        state.manualPhase = state.manualPhase(keepMask, keepMask);
        state.syntheticDistanceRatio = ...
            state.syntheticDistanceRatio(keepMask, keepMask);
        reconcileDisorderRealization();
        state.selectedSiteUids(ismember(state.selectedSiteUids, uids)) = [];
        if any(ismember(state.selectedEdgeUids, uids))
            state.selectedEdgeUids = [];
        end
        preserveCutoffDistance(previousCutoffDistance);
        recomputeAutomaticCouplings();
        refreshAll();
    end

    function siteDetuningChanged(src, ~)
        % Direct element edits convert a generated preset in place to Custom.
        if numel(state.selectedSiteUids) ~= 1
            return;
        end
        index = siteIndexFromUid(state.selectedSiteUids);
        if ~isfinite(src.Value)
            src.Value = state.nodes(index).detuning;
            return;
        end
        if src.Value == state.nodes(index).detuning
            return;
        end
        recordUndoState();
        convertPresetToCustomForEdit();
        state.nodes(index).detuning = src.Value;
        refreshAll();
    end

    function siteXChanged(src, ~)
        % Route X edits through the shared collision-aware position update.
        updateSelectedSitePosition(src.Value, true);
    end

    function siteYChanged(src, ~)
        % Route Y edits through the shared collision-aware position update.
        updateSelectedSitePosition(src.Value, false);
    end

    function updateSelectedSitePosition(value, updateX)
        % Preserve the absolute cutoff unless the new shortest edge exceeds it.
        if numel(state.selectedSiteUids) ~= 1
            return;
        end
        uid = state.selectedSiteUids;
        index = siteIndexFromUid(uid);
        if isempty(index)
            return;
        end
        proposed = [state.nodes(index).x, state.nodes(index).y];
        if updateX
            proposed(1) = value;
        else
            proposed(2) = value;
        end
        otherIndices = setdiff(1:numel(state.nodes), index);
        if ~isempty(otherIndices)
            otherCoordinates = [[state.nodes(otherIndices).x].', [state.nodes(otherIndices).y].'];
            if any(vecnorm(otherCoordinates - proposed, 2, 2) < 1e-8)
                refreshSelectionEditor();
                return;
            end
        end
        if isequaln(proposed, [state.nodes(index).x, state.nodes(index).y])
            return;
        end
        previousCutoffDistance = state.cutoffDistance;
        recordUndoState();
        convertPresetToCustomForEdit();
        state.nodes(index).x = proposed(1);
        state.nodes(index).y = proposed(2);
        reorderSites();
        preserveCutoffDistance(previousCutoffDistance);
        recomputeAutomaticCouplings();
        removeUnavailableEdgeSelection();
        refreshAll();
    end

    %% Directional coupling editing

    function directionalDistanceChanged(~, ~, isForward)
        % Distance is authoritative here: evaluate its ratio through the profile.
        [firstIndex, secondIndex] = selectedEdgeIndices();
        if isempty(firstIndex)
            return;
        end
        if isForward
            distanceField = ui.forwardDistanceField;
            sourceIndex = firstIndex;
            destinationIndex = secondIndex;
        else
            distanceField = ui.backwardDistanceField;
            sourceIndex = secondIndex;
            destinationIndex = firstIndex;
        end
        distance = distanceField.Value;
        dMin = minimumSiteDistance();
        if isnan(distance) || distance < 0 || isempty(dMin)
            refreshSelectionEditor();
            return;
        end
        currentDistance = displayedEffectiveDistance( ...
            sourceIndex, destinationIndex, dMin);
        if displayedScalarsMatch(distance, currentDistance)
            refreshSelectionEditor();
            return;
        end
        distanceRatio = distance/dMin;
        magnitude = geometricCoupling(distanceRatio);
        if ~isfinite(magnitude) || magnitude < 0
            refreshSelectionEditor();
            return;
        end
        phase = displayedPhase(sourceIndex, destinationIndex);
        previousValue = state.couplings(sourceIndex, destinationIndex);
        previousManual = state.manual(sourceIndex, destinationIndex);
        previousPhase = state.manualPhase(sourceIndex, destinationIndex);
        previousDistanceRatio = ...
            state.syntheticDistanceRatio(sourceIndex, destinationIndex);
        newValue = magnitude*exp(1i*phase);
        willChange = ~previousManual ...
            || ~isequaln(previousValue, newValue) ...
            || ~isequaln(previousPhase, phase) ...
            || ~isequaln(previousDistanceRatio, distanceRatio);
        if ~willChange
            refreshSelectionEditor();
            return;
        end
        recordUndoState();
        convertPresetToCustomForEdit();
        state.couplings(sourceIndex, destinationIndex) = newValue;
        state.manual(sourceIndex, destinationIndex) = true;
        state.manualPhase(sourceIndex, destinationIndex) = phase;
        state.syntheticDistanceRatio(sourceIndex, destinationIndex) = distanceRatio;
        refreshAll();
    end

    function directionalPolarChanged(~, ~, isForward, magnitudeChanged)
        % Magnitude edits update the profile-equivalent ratio; phase-only edits do not.
        [firstIndex, secondIndex] = selectedEdgeIndices();
        if isempty(firstIndex)
            return;
        end
        if isForward
            magnitudeField = ui.forwardMagnitudeField;
            phaseField = ui.forwardPhaseField;
            sourceIndex = firstIndex;
            destinationIndex = secondIndex;
        else
            magnitudeField = ui.backwardMagnitudeField;
            phaseField = ui.backwardPhaseField;
            sourceIndex = secondIndex;
            destinationIndex = firstIndex;
        end
        magnitude = magnitudeField.Value;
        phaseDegrees = phaseField.Value;
        phase = deg2rad(phaseDegrees);
        if ~isfinite(magnitude) || magnitude < 0 ...
                || ~isfinite(phaseDegrees) || ~isfinite(phase)
            refreshSelectionEditor();
            return;
        end
        wasManual = state.manual(sourceIndex, destinationIndex);
        previousValue = state.couplings(sourceIndex, destinationIndex);
        previousPhase = state.manualPhase(sourceIndex, destinationIndex);
        previousSyntheticDistanceRatio = ...
            state.syntheticDistanceRatio(sourceIndex, destinationIndex);
        currentPhaseDegrees = rad2deg(displayedPhase( ...
            sourceIndex, destinationIndex));
        if displayedScalarsMatch(magnitude, abs(previousValue)) ...
                && displayedScalarsMatch(phaseDegrees, currentPhaseDegrees)
            refreshSelectionEditor();
            return;
        end
        value = magnitude*exp(1i*phase);
        actualMagnitudeChange = magnitudeChanged ...
            && magnitude ~= abs(previousValue);
        proposedDistanceRatio = previousSyntheticDistanceRatio;
        if actualMagnitudeChange || ~wasManual || isnan(previousSyntheticDistanceRatio)
            proposedDistanceRatio = inferSyntheticDistanceRatio(magnitude);
            if ~couplingMatchesDistanceRatio(magnitude, proposedDistanceRatio)
                refreshSelectionEditor();
                return;
            end
        end
        willChange = ~wasManual ...
            || ~isequaln(previousValue, value) ...
            || ~isequaln(previousPhase, phase) ...
            || ~isequaln(previousSyntheticDistanceRatio, proposedDistanceRatio);
        if ~willChange
            refreshSelectionEditor();
            return;
        end
        recordUndoState();
        convertPresetToCustomForEdit();
        state.couplings(sourceIndex, destinationIndex) = value;
        state.manual(sourceIndex, destinationIndex) = true;
        state.manualPhase(sourceIndex, destinationIndex) = phase;
        state.syntheticDistanceRatio(sourceIndex, destinationIndex) = proposedDistanceRatio;
        refreshAll();
    end

    %% Geometric coupling model and site ordering

    function tf = couplingMatchesDistanceRatio(magnitude, distanceRatio)
        % Compare a manual magnitude with the active law at a stored ratio.
        reconstructedMagnitude = geometricCoupling(distanceRatio);
        if ~isfinite(reconstructedMagnitude) || reconstructedMagnitude < 0
            tf = false;
            return;
        end
        scale = max([magnitude, reconstructedMagnitude, realmin]);
        tf = abs(reconstructedMagnitude - magnitude) <= 1e-10*scale;
    end

    function tf = displayedScalarsMatch(firstValue, secondValue)
        % Ignore numerical roundoff when a displayed property is re-entered.
        if isequaln(firstValue, secondValue)
            tf = true;
            return;
        end
        if ~isfinite(firstValue) || ~isfinite(secondValue)
            tf = false;
            return;
        end
        scale = max([abs(firstValue), abs(secondValue), 1]);
        tf = abs(firstValue - secondValue) <= 1e-10*scale;
    end

    function restoreGeometricCouplings(~, ~)
        % An edge restores both directions; a node restores every incident direction.
        count = numel(state.nodes);
        restoreMask = false(count);
        if isscalar(state.selectedSiteUids)
            siteIndex = siteIndexFromUid(state.selectedSiteUids);
            if isempty(siteIndex)
                return;
            end
            restoreMask(siteIndex, :) = true;
            restoreMask(:, siteIndex) = true;
            restoreMask(siteIndex, siteIndex) = false;
        else
            [firstIndex, secondIndex] = selectedEdgeIndices();
            if isempty(firstIndex)
                return;
            end
            restoreMask(firstIndex, secondIndex) = true;
            restoreMask(secondIndex, firstIndex) = true;
        end
        manualRestoreMask = restoreMask & state.manual;
        if ~any(manualRestoreMask, 'all')
            return;
        end
        recordUndoState();
        convertPresetToCustomForEdit();
        state.manual(manualRestoreMask) = false;
        state.manualPhase(manualRestoreMask) = nan;
        state.syntheticDistanceRatio(manualRestoreMask) = nan;
        recomputeAutomaticCouplings();
        removeUnavailableEdgeSelection();
        refreshAll();
    end

    function distanceRatio = inferSyntheticDistanceRatio(magnitude)
        % Invert |kappa| to q=d_eff/d_min; zero coupling maps to infinite q.
        if magnitude == 0
            distanceRatio = Inf;
        else
            switch state.couplingProfile
                case "Power law"
                    distanceRatio = exp(-log(magnitude)/state.gamma);
                case "Exponential"
                    distanceRatio = max(0, 1 - log(magnitude)/state.gamma);
                case "Fiber asymptotic"
                    distanceRatio = inverseFiberDistanceRatio(magnitude);
                otherwise
                    error('lattice_config:UnknownCouplingProfile', ...
                        'Unknown coupling profile "%s".', state.couplingProfile);
            end
        end
    end

    function distanceRatio = inverseFiberDistanceRatio(magnitude)
        % Solve gamma*(q-1)+log(q)/2+log(|kappa|)=0 in log(q) by bisection.
        logMagnitude = log(magnitude);
        if logMagnitude == 0
            distanceRatio = 1;
            return;
        elseif logMagnitude > 0
            lowerLogRatio = -2*(logMagnitude + 1);
            upperLogRatio = 0;
        else
            lowerLogRatio = 0;
            upperLogRatio = 2*(-logMagnitude + 1);
        end
        for iteration = 1:80
            middleLogRatio = (lowerLogRatio + upperLogRatio)/2;
            residual = state.gamma*expm1(middleLogRatio) ...
                + 0.5*middleLogRatio + logMagnitude;
            if residual < 0
                lowerLogRatio = middleLogRatio;
            else
                upperLogRatio = middleLogRatio;
            end
        end
        distanceRatio = exp((lowerLogRatio + upperLogRatio)/2);
    end

    function reorderSites()
        % Group nearly equal Y values into rows and permute all matrices together.
        count = numel(state.nodes);
        if count < 2
            return;
        end
        dMin = minimumSiteDistance();
        rowTolerance = 1e-10;
        if ~isempty(dMin)
            rowTolerance = rowTolerance*max(1, dMin);
        end
        rowKey = round([state.nodes.y].'/rowTolerance)*rowTolerance;
        orderingData = [-rowKey, [state.nodes.x].', [state.nodes.uid].'];
        [~, order] = sortrows(orderingData, [1 2 3]);
        state.nodes = state.nodes(order);
        state.couplings = state.couplings(order, order);
        state.manual = state.manual(order, order);
        state.manualPhase = state.manualPhase(order, order);
        state.syntheticDistanceRatio = state.syntheticDistanceRatio(order, order);
    end

    function recomputeAutomaticCouplings()
        % Rebuild only unmasked directions; manual overrides bypass the cutoff.
        count = numel(state.nodes);
        if count == 0
            state.couplings = complex(zeros(0));
            return;
        end
        state.couplings(1:count + 1:end) = 0;
        dMin = minimumSiteDistance();
        if isempty(dMin)
            state.couplings(~state.manual) = 0;
            return;
        end
        cutoffDistance = state.cutoffDistance;
        for firstIndex = 1:count - 1
            firstPosition = [state.nodes(firstIndex).x, state.nodes(firstIndex).y];
            for secondIndex = firstIndex + 1:count
                secondPosition = [state.nodes(secondIndex).x, state.nodes(secondIndex).y];
                separation = norm(secondPosition - firstPosition);
                distanceRatio = separation/dMin;
                if abs(distanceRatio - 1) <= 1e-10
                    distanceRatio = 1;
                end
                geometricValue = geometricCoupling(distanceRatio);
                if separation > cutoffDistance*(1 + 1e-10)
                    geometricValue = 0;
                end
                if ~state.manual(firstIndex, secondIndex)
                    state.couplings(firstIndex, secondIndex) = geometricValue;
                end
                if ~state.manual(secondIndex, firstIndex)
                    state.couplings(secondIndex, firstIndex) = geometricValue;
                end
            end
        end
    end

    function value = geometricCoupling(distanceRatio)
        % All profiles are normalized to one at distanceRatio = 1.
        switch state.couplingProfile
            case "Power law"
                value = distanceRatio^(-state.gamma);
            case "Exponential"
                value = exp(-state.gamma*(distanceRatio - 1));
            case "Fiber asymptotic"
                value = distanceRatio^(-0.5) ...
                    * exp(-state.gamma*(distanceRatio - 1));
            otherwise
                error('lattice_config:UnknownCouplingProfile', ...
                    'Unknown coupling profile "%s".', state.couplingProfile);
        end
    end

    function cutoffDistance = correctedCutoffDistance(cutoffDistance)
        % Enforce a valid absolute radius for the current geometry.
        dMin = minimumSiteDistance();
        if ~isempty(dMin) && cutoffDistance < dMin
            cutoffDistance = dMin + eps(dMin);
        end
    end

    function preserveCutoffDistance(cutoffDistance)
        % Restore a drag/edit baseline, raising it only for a larger d_min.
        if isempty(cutoffDistance)
            return;
        end
        state.cutoffDistance = correctedCutoffDistance(cutoffDistance);
    end

    function dMin = minimumSiteDistance()
        % Return the shortest nonzero pair separation, or [] below two sites.
        count = numel(state.nodes);
        if count < 2
            dMin = [];
            return;
        end
        dMin = Inf;
        for firstIndex = 1:count - 1
            firstPosition = [state.nodes(firstIndex).x, state.nodes(firstIndex).y];
            for secondIndex = firstIndex + 1:count
                secondPosition = [state.nodes(secondIndex).x, state.nodes(secondIndex).y];
                separation = norm(secondPosition - firstPosition);
                if separation > 1e-12
                    dMin = min(dMin, separation);
                end
            end
        end
        if ~isfinite(dMin)
            dMin = [];
        end
    end

    %% Interface refresh and control layout

    function refreshAll()
        % Synchronize controls, graph, and published state from one snapshot.
        refreshCouplingEquation();
        refreshPresetControls();
        refreshSelectionEditor();
        redrawGraph();
    end

    function refreshCouplingEquation()
        % Keep the displayed formula identical to the evaluated profile.
        ui.profileDropdown.Value = char(state.couplingProfile);
        ui.formulaLabel.Text = couplingFormula(state.couplingProfile);
    end

    function formula = couplingFormula(profile)
        % Return the compact LaTeX form used by the sublattice property panel.
        switch string(profile)
            case "Power law"
                formula = ...
                    ['$\kappa_{i,j}=(1+r_{i,j})^{-\gamma}' ...
                     '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
            case "Exponential"
                formula = ...
                    ['$\kappa_{i,j}=\exp(-\gamma r_{i,j})' ...
                     '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
            case "Fiber asymptotic"
                formula = ...
                    ['$\kappa_{i,j}=(1+r_{i,j})^{-1/2}\exp(-\gamma r_{i,j})' ...
                     '\:\:\:,\quad r_{i,j}=d_{i,j}/d_{\min}-1$'];
            otherwise
                error('lattice_config:UnknownCouplingProfile', ...
                    'Unknown coupling profile "%s".', string(profile));
        end
    end

    function refreshPresetControls()
        % Refresh common controls and the active family's compact fields.
        isCustom = state.preset == "Custom";
        ui.presetDropdown.Value = char(state.preset);
        ui.gammaField.Value = state.gamma;
        ui.cutoffField.Value = state.cutoffDistance;
        ui.disorderDropdown.Value = char(state.disorder.type);
        ui.disorderStrengthField.Value = state.disorder.strength;
        refreshDisorderStrengthControl( ...
            ui.disorderStrengthField, state.disorder.type);
        leftGrid.RowHeight = num2cell(leftRowHeights);
        dMin = minimumSiteDistance();
        if ~isempty(dMin)
            cutoffLabel.Text = sprintf('Range cutoff (dₘᵢₙ=%.3f)', dMin);
        else
            cutoffLabel.Text = 'Range cutoff';
        end
        if isCustom
            setVisibility(ui.presetParameterCells, false);
        else
            refreshPresetParameterControls();
        end
        ui.acceptButton.Enable = 'on';
    end

    function refreshDisorderStrengthControl(field, disorderType)
        % Disable and gray the strength field when disorder is inactive.
        enabled = string(disorderType) ~= "No disorder";
        if enabled
            field.Enable = 'on';
            field.Editable = 'on';
            field.BackgroundColor = [1 1 1];
            field.FontColor = [0 0 0];
        else
            field.Value = 0;
            field.Editable = 'off';
            field.Enable = 'off';
            field.BackgroundColor = readOnlyBackgroundColor;
            field.FontColor = readOnlyFontColor;
        end
    end

    function refreshPresetParameterControls()
        % Fixed field anchors make labels follow boxes rather than shift them.
        setVisibility(ui.presetParameterCells, false);
        [keys, labels, tags, integerMask, lowerBounds] = ...
            presetParameterDefinition(state.preset);
        parameterRows = presetParameterLayout(state.preset);
        settings = state.presetSettings.(char(state.preset));
        usesMiddleRow = any(state.preset == ["Triangular", "Ring"]);
        % Sparse two-parameter families share one row centered across rows 1–2.
        presetParameterMiddleGrid.Visible = ...
            matlab.lang.OnOffSwitchState(usesMiddleRow);
        for parameterRow = 1:2
            presetParameterRowGrids(parameterRow).Visible = ...
                matlab.lang.OnOffSwitchState(~usesMiddleRow);
        end
        if usesMiddleRow
            middleSlots = [1 3];
            for parameterIndex = 1:numel(keys)
                slot = middleSlots(parameterIndex);
                ui.presetParameterCells(parameterIndex).Parent = ...
                    presetParameterMiddleGrid;
                ui.presetParameterCells(parameterIndex).Layout.Row = 2;
                ui.presetParameterCells(parameterIndex).Layout.Column = ...
                    [2*slot - 1, 2*slot];
            end
        else
            for parameterRow = 1:2
                rowParameterIndices = find(parameterRows == parameterRow);
                [rowColumnWidths, rowSlots] = presetParameterRowLayout( ...
                    state.preset, parameterRow, numel(rowParameterIndices));
                presetParameterRowGrids(parameterRow).ColumnWidth = ...
                    rowColumnWidths;
                for rowPosition = 1:numel(rowParameterIndices)
                    parameterIndex = rowParameterIndices(rowPosition);
                    slot = rowSlots(rowPosition);
                    ui.presetParameterCells(parameterIndex).Parent = ...
                        presetParameterRowGrids(parameterRow);
                    ui.presetParameterCells(parameterIndex).Layout.Row = 1;
                    ui.presetParameterCells(parameterIndex).Layout.Column = ...
                        [2*slot - 1, 2*slot];
                end
            end
        end
        for parameterIndex = 1:numel(keys)
            ui.presetParameterLabels(parameterIndex).Text = labels{parameterIndex};
            ui.presetParameterFields(parameterIndex).Tag = tags{parameterIndex};
            ui.presetParameterFields(parameterIndex).Limits = [-Inf Inf];
            ui.presetParameterFields(parameterIndex).RoundFractionalValues = 'off';
            ui.presetParameterFields(parameterIndex).Value = ...
                settings.(keys{parameterIndex});
            ui.presetParameterFields(parameterIndex).Limits = ...
                [lowerBounds(parameterIndex) Inf];
            ui.presetParameterFields(parameterIndex).RoundFractionalValues = ...
                matlab.lang.OnOffSwitchState(integerMask(parameterIndex));
            ui.presetParameterCells(parameterIndex).Visible = 'on';
        end
    end

    function [columnWidths, slots] = presetParameterRowLayout( ...
            preset, parameterRow, parameterCount)
        % Return field-first grid weights for the active family and row.
        usesThirdAnchors = preset == "Lieb";
        usesInsetOuterAnchors = ...
            any(preset == ["Square", "Jx"]) && parameterRow == 1;
        centersSingleField = any(preset == ["Square", "Jx"]) && ...
            parameterRow == 2;
        if usesInsetOuterAnchors
            columnWidths = {'7x', presetParameterFieldWidth, '2x', ...
                presetParameterFieldWidth, '2x', ...
                presetParameterFieldWidth, '7x'};
            slots = [1 3];
            return
        end
        if usesThirdAnchors || centersSingleField
            columnWidths = {'1x', presetParameterFieldWidth, '2x', ...
                presetParameterFieldWidth, '2x', ...
                presetParameterFieldWidth, '1x'};
            if centersSingleField
                slots = 2;
            elseif parameterCount == 2
                slots = [1 3];
            else
                slots = 1:parameterCount;
            end
            return
        end
        columnWidths = repmat({0}, 1, 7);
        for slot = 1:parameterCount
            columnWidths{2*slot - 1} = '1x';
            columnWidths{2*slot} = presetParameterFieldWidth;
        end
        columnWidths{2*parameterCount + 1} = '1x';
        slots = 1:parameterCount;
    end

    function refreshSelectionEditor()
        % A preset selection temporarily replaces family controls with element fields.
        hideSelectionControls();
        isCustom = state.preset == "Custom";
        hasElementSelection = ~isempty(state.selectedSiteUids) ...
            || numel(state.selectedEdgeUids) == 2;
        showElementProperties = isCustom || hasElementSelection;
        selectionGrid.Visible = ...
            matlab.lang.OnOffSwitchState(showElementProperties);
        presetPropertiesGrid.Visible = ...
            matlab.lang.OnOffSwitchState(~showElementProperties);
        if ~showElementProperties
            ui.selectionTitleLabel.Text = ...
                sprintf('%s properties', char(state.preset));
            return;
        end
        % Element fields always expose the disorder-free base model. Merely
        % viewing them does not convert a generated preset to Custom.
        ui.selectionTitleLabel.Text = 'Element properties';
        selectionGrid.RowHeight = {'1x', 0, 0, 0};
        if isscalar(state.selectedSiteUids)
            siteIndex = siteIndexFromUid(state.selectedSiteUids);
            ui.selectionTitleLabel.Text = sprintf('Node %d properties', siteIndex);
            selectionGrid.RowHeight = {26, 26, 26, 26};
            ui.detuningField.Value = state.nodes(siteIndex).detuning;
            ui.xField.Value = state.nodes(siteIndex).x;
            ui.yField.Value = state.nodes(siteIndex).y;
            setVisibility([ui.detuningLabel, ui.detuningField, ui.xLabel, ui.xField, ...
                ui.yLabel, ui.yField, ui.restoreButton], true);
            return;
        elseif numel(state.selectedSiteUids) > 1
            ui.selectionTitleLabel.Text = sprintf( ...
                '%d nodes selected', numel(state.selectedSiteUids));
            ui.selectionHint.Text = 'Drag any selected site to move the group';
            ui.selectionHint.Visible = 'on';
            return;
        end

        [firstIndex, secondIndex] = selectedEdgeIndices();
        if ~isempty(firstIndex)
            ui.selectionTitleLabel.Text = sprintf( ...
                'Edge %d↔%d properties', firstIndex, secondIndex);
            selectionGrid.RowHeight = {26, 26, 26, 26};
            ui.forwardLabel.Text = sprintf('%d → %d', firstIndex, secondIndex);
            ui.backwardLabel.Text = sprintf('%d ← %d', firstIndex, secondIndex);
            forwardValue = state.couplings(firstIndex, secondIndex);
            backwardValue = state.couplings(secondIndex, firstIndex);
            dMin = minimumSiteDistance();
            ui.forwardDistanceField.Value = displayedEffectiveDistance( ...
                firstIndex, secondIndex, dMin);
            ui.backwardDistanceField.Value = displayedEffectiveDistance( ...
                secondIndex, firstIndex, dMin);
            ui.forwardMagnitudeField.Value = abs(forwardValue);
            ui.backwardMagnitudeField.Value = abs(backwardValue);
            ui.forwardPhaseField.Value = rad2deg(displayedPhase(firstIndex, secondIndex));
            ui.backwardPhaseField.Value = rad2deg(displayedPhase(secondIndex, firstIndex));
            setVisibility([ui.directionHeader, ui.distanceHeader, ui.magnitudeHeader, ...
                ui.phaseHeader, ui.forwardLabel, ui.forwardDistanceField, ...
                ui.forwardMagnitudeField, ui.forwardPhaseField, ui.backwardLabel, ...
                ui.backwardDistanceField, ui.backwardMagnitudeField, ...
                ui.backwardPhaseField, ui.restoreButton], true);
            return;
        end

        ui.selectionHint.Text = 'Select a site or edge in the editing plane';
        ui.selectionHint.Visible = 'on';
    end

    function phase = displayedPhase(sourceIndex, destinationIndex)
        % Prefer a stored manual phase; otherwise infer it from the base value.
        phase = state.manualPhase(sourceIndex, destinationIndex);
        if isnan(phase)
            phase = angle(state.couplings(sourceIndex, destinationIndex));
        end
    end

    function distance = displayedEffectiveDistance(sourceIndex, destinationIndex, dMin)
        % Convert the direction's synthetic ratio into an editable distance.
        distanceRatio = state.syntheticDistanceRatio(sourceIndex, destinationIndex);
        if isnan(distanceRatio)
            distanceRatio = inferSyntheticDistanceRatio( ...
                abs(state.couplings(sourceIndex, destinationIndex)));
        end
        if isempty(dMin)
            distance = nan;
        else
            distance = dMin*distanceRatio;
        end
    end

    function hideSelectionControls()
        % Reset every mutually exclusive node/edge property control to hidden.
        ui.selectionHint.Visible = 'off';
        setVisibility([ui.detuningLabel, ui.detuningField, ui.xLabel, ui.xField, ...
            ui.yLabel, ui.yField, ui.directionHeader, ui.distanceHeader, ...
            ui.magnitudeHeader, ui.phaseHeader, ui.forwardLabel, ...
            ui.forwardDistanceField, ui.forwardMagnitudeField, ui.forwardPhaseField, ...
            ui.backwardLabel, ui.backwardDistanceField, ui.backwardMagnitudeField, ...
            ui.backwardPhaseField, ui.restoreButton], false);
    end

    function setVisibility(handles, isVisible)
        % Apply one visibility state to a heterogeneous graphics-handle list.
        value = matlab.lang.OnOffSwitchState(isVisible);
        for handleIndex = 1:numel(handles)
            handles(handleIndex).Visible = value;
        end
    end

    %% Graph rendering and edge geometry

    function redrawGraph()
        % Render effective disorder values while preserving base values for editing.
        if ~isvalid(ui.axes)
            return;
        end
        delete(ui.axes.Children);
        hold(ui.axes, 'on');

        count = numel(state.nodes);
        if count == 0
            text(ui.axes, 0, 0, 'Double-click to place the first site', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', regularFontSize + 4, ...
                'Color', [0.55 0.58 0.62], ...
                'HitTest', 'off');
            hold(ui.axes, 'off');
            return;
        end

        dMin = minimumSiteDistance();
        effectiveCouplingValues = effectiveCouplings();
        effectiveDetuningValues = effectiveDetunings();
        for firstIndex = 1:count - 1
            for secondIndex = firstIndex + 1:count
                firstBaseValue = state.couplings(firstIndex, secondIndex);
                secondBaseValue = state.couplings(secondIndex, firstIndex);
                firstValue = effectiveCouplingValues(firstIndex, secondIndex);
                secondValue = effectiveCouplingValues(secondIndex, firstIndex);
                firstManual = state.manual(firstIndex, secondIndex);
                secondManual = state.manual(secondIndex, firstIndex);
                if abs(firstBaseValue) <= 1e-14 && abs(secondBaseValue) <= 1e-14 ...
                        && ~firstManual && ~secondManual
                    continue;
                end
                uidA = state.nodes(firstIndex).uid;
                uidB = state.nodes(secondIndex).uid;
                isSelected = selectedEdgeMatches(uidA, uidB);
                comparisonScale = max(abs([firstValue, secondValue]));
                valuesMatch = comparisonScale == 0 || ...
                    abs(firstValue - secondValue) <= 1e-12*comparisonScale;
                if valuesMatch
                    % Reciprocal equal directions share one centered path.
                    if state.disorder.type == "Coupling disorder"
                        syntheticDistanceRatio = mean([ ...
                            displayedSyntheticDistanceRatio( ...
                                firstIndex, secondIndex, firstValue), ...
                            displayedSyntheticDistanceRatio( ...
                                secondIndex, firstIndex, secondValue)]);
                    else
                        syntheticDistanceRatio = nan;
                        manualRatios = [ ...
                            state.syntheticDistanceRatio(firstIndex, secondIndex), ...
                            state.syntheticDistanceRatio(secondIndex, firstIndex)];
                        manualRatios = manualRatios(~isnan(manualRatios));
                        if ~isempty(manualRatios)
                            syntheticDistanceRatio = mean(manualRatios);
                        end
                    end
                    drawEdgeDirection(firstIndex, secondIndex, firstValue, [0 0 0], ...
                        scaledSyntheticDistance(syntheticDistanceRatio, dMin), ...
                        0, false, uidA, uidB, isSelected);
                else
                    % Offset asymmetric paths, then draw both arrows on top.
                    [firstColor, secondColor] = ...
                        couplingAsymmetryColors(firstValue, secondValue);
                    [firstArrowPath, firstArrowAlpha] = ...
                        drawEdgeDirection(firstIndex, secondIndex, firstValue, firstColor, ...
                        scaledSyntheticDistance( ...
                            displayedSyntheticDistanceRatio( ...
                                firstIndex, secondIndex, firstValue), dMin), ...
                        0.13, true, ...
                        uidA, uidB, isSelected);
                    [secondArrowPath, secondArrowAlpha] = ...
                        drawEdgeDirection(secondIndex, firstIndex, secondValue, secondColor, ...
                        scaledSyntheticDistance( ...
                            displayedSyntheticDistanceRatio( ...
                                secondIndex, firstIndex, secondValue), dMin), ...
                        0.13, true, ...
                        uidA, uidB, isSelected);
                    callback = @(src, event) edgeClicked(src, event, uidA, uidB);
                    if ~isempty(firstArrowPath)
                        drawArrowhead(firstArrowPath, firstColor, firstArrowAlpha, ...
                            callback);
                    end
                    if ~isempty(secondArrowPath)
                        drawArrowhead(secondArrowPath, secondColor, secondArrowAlpha, ...
                            callback);
                    end
                end
            end
        end

        for siteIndex = 1:count
            uid = state.nodes(siteIndex).uid;
            position = [state.nodes(siteIndex).x, state.nodes(siteIndex).y];
            isSelected = any(state.selectedSiteUids == uid);
            detuning = effectiveDetuningValues(siteIndex);
            faceColor = detuningNodeColor(detuning);
            if isSelected
                selectionHalo = scatter(ui.axes, position(1), position(2), 405, ...
                    'o', 'filled', ...
                    'MarkerFaceColor', [0.98 0.70 0.25], ...
                    'MarkerEdgeColor', [0.62 0.32 0.04], ...
                    'LineWidth', 1.4, ...
                    'Tag', 'SublatticeSiteSelectionHalo', ...
                    'HitTest', 'off');
                selectionHalo.PickableParts = 'none';
            end
            marker = scatter(ui.axes, position(1), position(2), 310, ...
                'o', 'filled', ...
                'MarkerFaceColor', faceColor, ...
                'MarkerEdgeColor', [0 0 0], ...
                'LineWidth', 1.8, ...
                'Tag', 'SublatticeSite', ...
                'ButtonDownFcn', @(src, event) nodeClicked(src, event, uid));
            marker.PickableParts = 'all';
            text(ui.axes, position(1), position(2), sprintf('%d', siteIndex), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'FontSize', regularFontSize - 2, ...
                'Color', [1 1 1], ...
                'HitTest', 'on', ...
                'PickableParts', 'all', ...
                'ButtonDownFcn', @(src, event) nodeClicked(src, event, uid));
        end
        drawSelectionRectangle();
        hold(ui.axes, 'off');
    end

    function distance = scaledSyntheticDistance(distanceRatio, dMin)
        % Convert a stored ratio to the path length represented in the graph.
        if isnan(distanceRatio) || isempty(dMin)
            distance = nan;
        else
            distance = dMin*distanceRatio;
        end
    end

    function drawSelectionRectangle()
        % Overlay the active box-selection marquee after graph redraws.
        if state.interaction.mode ~= "box" || ~state.interaction.moved
            return;
        end
        xBounds = sort([state.interaction.startPoint(1), state.interaction.currentPoint(1)]);
        yBounds = sort([state.interaction.startPoint(2), state.interaction.currentPoint(2)]);
        verticesX = [xBounds(1), xBounds(2), xBounds(2), xBounds(1)];
        verticesY = [yBounds(1), yBounds(1), yBounds(2), yBounds(2)];
        selectionPatch = patch(ui.axes, verticesX, verticesY, [0.20 0.48 0.78], ...
            'FaceAlpha', 0.08, ...
            'EdgeColor', [0.20 0.48 0.78], ...
            'LineStyle', '--', ...
            'LineWidth', 1.2, ...
            'Tag', 'SublatticeSelectionBox', ...
            'HitTest', 'off', ...
            'PickableParts', 'none');
        selectionPatch.Clipping = 'on';
    end

    function [arrowPath, arrowAlpha] = drawEdgeDirection( ...
            sourceIndex, destinationIndex, value, color, syntheticDistance, ...
            offset, directed, uidA, uidB, isSelected)
        % Longer effective distance curves the path; shorter distance dots its core.
        arrowPath = zeros(0, 2);
        arrowAlpha = 0;
        source = [state.nodes(sourceIndex).x, state.nodes(sourceIndex).y];
        destination = [state.nodes(destinationIndex).x, state.nodes(destinationIndex).y];
        delta = destination - source;
        actualDistance = norm(delta);
        if actualDistance <= 1e-12
            return;
        end
        magnitude = abs(value);
        direction = delta/actualDistance;
        normal = [-direction(2), direction(1)];
        source = source + offset*normal;
        destination = destination + offset*normal;

        callback = @(src, event) edgeClicked(src, event, uidA, uidB);
        if magnitude <= 1e-14
            path = [source; destination];
            drawPath(path, [0.55 0.58 0.62], 0, 9, '-', true, callback);
            if isSelected
                drawPath(path, [0.98 0.67 0.12], 0.12, 5.5, ':', false, callback);
            else
                drawPath(path, color, 0.08, 1.2, ':', false, callback);
            end
            if directed
                arrowPath = path;
                arrowAlpha = 0.08;
            end
            return;
        end

        if isnan(syntheticDistance)
            effectiveDistance = actualDistance;
        elseif isinf(syntheticDistance)
            effectiveDistance = min(4*actualDistance, 24);
        else
            effectiveDistance = syntheticDistance;
        end

        tolerance = 1e-6*max(1, actualDistance);
        isCurved = effectiveDistance > actualDistance + tolerance;
        plotBox = axesPlotBoxPosition();
        pixelsPerUnit = min(plotBox(3)/diff(ui.axes.XLim), ...
            plotBox(4)/diff(ui.axes.YLim));
        minimumVisibleShortening = max(0.02*actualDistance, ...
            6/max(pixelsPerUnit, eps));
        % Suppress sub-pixel shortcut ticks that look like stray crosses.
        isShortened = actualDistance - effectiveDistance ...
            > minimumVisibleShortening;
        if isCurved
            path = curvedPath(source, destination, effectiveDistance);
        else
            path = [source; destination];
        end

        if isSelected
            drawPath(path, [0.98 0.67 0.12], 0.16, 6.5, '-', false, callback);
        end

        alpha = 1 - exp(-2*magnitude);
        width = 1.3 + 1.4*min(1, sqrt(magnitude));

        drawPath(path, color, 0, max(9, width + 6), '-', true, callback);
        if isShortened
            solidFraction = max(0, min(1, effectiveDistance/actualDistance));
            firstBreak = solidFraction/2;
            secondBreak = 1 - firstBreak;
            drawPath(pathFraction(path, 0, firstBreak), ...
                color, alpha, width, '-', false, callback);
            drawPath(pathFraction(path, firstBreak, secondBreak), ...
                color, alpha, width, ':', false, callback);
            drawPath(pathFraction(path, secondBreak, 1), ...
                color, alpha, width, '-', false, callback);
            drawShortcutBoundary(path, firstBreak, color, alpha, width);
            drawShortcutBoundary(path, secondBreak, color, alpha, width);
        else
            drawPath(path, color, alpha, width, '-', false, callback);
        end

        if directed && magnitude > 0
            arrowPath = path;
            arrowAlpha = alpha;
        end

    end

    %% Drawing and path primitives

    function drawPath(path, color, alpha, width, style, isHitTarget, callback)
        % Draw either a visible pre-blended stroke or its transparent hit ribbon.
        if size(path, 1) < 2
            return;
        end
        if isHitTarget
            drawHitRibbon(path, callback);
            return;
        end
        background = ui.axes.Color;
        % Pre-blend against the axes because line objects do not expose alpha.
        displayColor = background + alpha*(color - background);
        line(ui.axes, path(:, 1), path(:, 2), ...
            'Color', displayColor, ...
            'LineWidth', width, ...
            'LineStyle', style, ...
            'Clipping', 'on', ...
            'Tag', 'SublatticeEdge', ...
            'HitTest', 'off', ...
            'PickableParts', 'none');
    end

    function drawHitRibbon(path, callback)
        % A transparent ribbon keeps thin and curved edges easy to select.
        tangents = zeros(size(path));
        tangents(1, :) = path(2, :) - path(1, :);
        tangents(end, :) = path(end, :) - path(end - 1, :);
        if size(path, 1) > 2
            tangents(2:end - 1, :) = path(3:end, :) - path(1:end - 2, :);
        end
        tangentNorms = vecnorm(tangents, 2, 2);
        tangentNorms(tangentNorms <= eps) = 1;
        tangents = tangents./tangentNorms;
        normals = [-tangents(:, 2), tangents(:, 1)];
        dataSpan = min(diff(ui.axes.XLim), diff(ui.axes.YLim));
        halfWidth = max(0.02, 0.0065*dataSpan);
        upperEdge = path + halfWidth*normals;
        lowerEdge = flipud(path - halfWidth*normals);
        vertices = [upperEdge; lowerEdge];
        graphic = patch(ui.axes, ...
            'XData', vertices(:, 1), ...
            'YData', vertices(:, 2), ...
            'FaceColor', [0 0 0], ...
            'FaceAlpha', 0.001, ...
            'EdgeColor', 'none', ...
            'Clipping', 'on', ...
            'Tag', 'SublatticeEdgeHitTarget', ...
            'ButtonDownFcn', callback);
        graphic.PickableParts = 'all';
        graphic.HitTest = 'on';
    end

    function drawShortcutBoundary(path, fraction, color, alpha, width)
        % Mark a solid-to-dotted transition with a short transverse tick.
        [point, tangent] = pathPointAndTangent(path, fraction);
        if isempty(point)
            return;
        end
        normal = [-tangent(2), tangent(1)];
        halfLength = 0.0065*min(diff(ui.axes.XLim), diff(ui.axes.YLim));
        endpoints = [point - halfLength*normal; point + halfLength*normal];
        background = ui.axes.Color;
        displayColor = background + alpha*(color - background);
        line(ui.axes, endpoints(:, 1), endpoints(:, 2), ...
            'Color', displayColor, ...
            'LineWidth', max(1.5, width), ...
            'LineStyle', '-', ...
            'Clipping', 'on', ...
            'Tag', 'SublatticeEdgeShortcutBoundary', ...
            'HitTest', 'off', ...
            'PickableParts', 'none');
    end

    function drawArrowhead(path, color, alpha, callback)
        % Place one fixed-pixel arrowhead at the path midpoint.
        if size(path, 1) < 2
            return;
        end
        [center, tangent] = pathPointAndTangent(path, 0.5);
        if isempty(center)
            return;
        end
        normal = [-tangent(2), tangent(1)];
        % Keep arrowheads compact in screen space regardless of zoom or domain.
        plotBox = axesPlotBoxPosition();
        pixelsPerUnit = min(plotBox(3)/diff(ui.axes.XLim), ...
            plotBox(4)/diff(ui.axes.YLim));
        tipLength = 8.5/max(pixelsPerUnit, eps);
        baseLength = 7.5/max(pixelsPerUnit, eps);
        halfWidth = 5/max(pixelsPerUnit, eps);
        tip = center + tipLength*tangent;
        base = center - baseLength*tangent;
        vertices = [tip; base + halfWidth*normal; base - halfWidth*normal];
        arrow = patch(ui.axes, ...
            'XData', vertices(:, 1), ...
            'YData', vertices(:, 2), ...
            'FaceColor', color, ...
            'FaceAlpha', max(0.65, alpha), ...
            'EdgeColor', 'none', ...
            'Tag', 'SublatticeEdgeArrow', ...
            'ButtonDownFcn', callback);
        arrow.PickableParts = 'all';
    end

    function [point, tangent] = pathPointAndTangent(path, fraction)
        % Interpolate by arc-length fraction so arrows and ticks stay centered.
        point = [];
        tangent = [];
        if size(path, 1) < 2
            return;
        end
        cumulative = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
        total = cumulative(end);
        if total <= eps
            return;
        end
        target = max(0, min(1, fraction))*total;
        if target <= 0
            lower = 1;
            upper = 2;
        elseif target >= total
            lower = size(path, 1) - 1;
            upper = size(path, 1);
        else
            upper = find(cumulative >= target, 1, 'first');
            lower = upper - 1;
        end
        segmentLength = cumulative(upper) - cumulative(lower);
        localFraction = (target - cumulative(lower))/max(segmentLength, eps);
        point = path(lower, :) + localFraction*(path(upper, :) - path(lower, :));
        tangent = path(upper, :) - path(lower, :);
        tangent = tangent/max(norm(tangent), eps);
    end

    function path = curvedPath(source, destination, targetLength)
        % Solve the sag of a quadratic arc whose length matches targetLength.
        chord = norm(destination - source);
        targetLength = max(chord, min(targetLength, max(4*chord, 24)));
        direction = (destination - source)/chord;
        normal = [-direction(2), direction(1)];
        low = 0;
        high = max(chord, 0.5);
        while quadraticArcLength(chord, high) < targetLength && high < 100
            high = 2*high;
        end
        for iteration = 1:45
            midpoint = (low + high)/2;
            if quadraticArcLength(chord, midpoint) < targetLength
                low = midpoint;
            else
                high = midpoint;
            end
        end
        sag = (low + high)/2;
        parameter = linspace(0, 1, 81).';
        straightPart = (1 - parameter).*source + parameter.*destination;
        path = straightPart + (4*sag*parameter.*(1 - parameter)).*normal;
    end

    function lengthValue = quadraticArcLength(chord, sag)
        % Evaluate the exact length of the quadratic arc used by curvedPath.
        if sag <= eps*max(1, chord)
            lengthValue = chord;
            return;
        end
        lengthValue = 0.5*sqrt(chord^2 + 16*sag^2) ...
            + chord^2/(8*sag)*asinh(4*sag/chord);
    end

    function segment = pathFraction(path, startFraction, endFraction)
        % Extract a path interval by arc length, not by vertex index.
        startFraction = max(0, min(1, startFraction));
        endFraction = max(startFraction, min(1, endFraction));
        cumulative = [0; cumsum(vecnorm(diff(path, 1, 1), 2, 2))];
        total = cumulative(end);
        if total <= eps || endFraction - startFraction <= eps
            segment = zeros(0, 2);
            return;
        end
        startDistance = startFraction*total;
        endDistance = endFraction*total;
        interior = cumulative > startDistance & cumulative < endDistance;
        startPoint = interpolatePath(path, cumulative, startDistance);
        endPoint = interpolatePath(path, cumulative, endDistance);
        segment = [startPoint; path(interior, :); endPoint];
    end

    function point = interpolatePath(path, cumulative, distance)
        % Interpolate a polyline at an absolute cumulative arc length.
        if distance <= 0
            point = path(1, :);
            return;
        elseif distance >= cumulative(end)
            point = path(end, :);
            return;
        end
        upper = find(cumulative >= distance, 1, 'first');
        lower = upper - 1;
        fraction = (distance - cumulative(lower))/ ...
            max(cumulative(upper) - cumulative(lower), eps);
        point = path(lower, :) + fraction*(path(upper, :) - path(lower, :));
    end

    %% Color mappings

    function [firstColor, secondColor] = couplingAsymmetryColors(firstValue, secondValue)
        % Color depends only on normalized directional magnitude imbalance.
        firstMagnitude = abs(firstValue);
        secondMagnitude = abs(secondValue);
        magnitudeSum = firstMagnitude + secondMagnitude;
        if magnitudeSum <= eps
            imbalance = 0;
        else
            imbalance = (firstMagnitude - secondMagnitude)/magnitudeSum;
        end
        firstColor = couplingImbalanceColor(imbalance);
        secondColor = couplingImbalanceColor(-imbalance);
    end

    function color = couplingImbalanceColor(imbalance)
        % Sharpen normalized imbalance into the black-to-green/red palette.
        neutral = [0 0 0];
        green = [0.10 0.72 0.27];
        red = [0.86 0.20 0.16];
        imbalance = max(-1, min(1, imbalance));
        if imbalance >= 0
            target = green;
        else
            target = red;
        end
        blend = min(1, 1.3*abs(imbalance));
        color = neutral + blend*(target - neutral);
    end

    function color = detuningNodeColor(detuning)
        % Saturate smoothly from black toward blue (positive) or orange (negative).
        neutral = [0.08 0.08 0.08];
        positive = [0.12 0.42 0.82];
        negative = [0.86 0.32 0.06];
        if detuning >= 0
            target = positive;
        else
            target = negative;
        end
        blend = 1 - exp(-abs(detuning));
        color = neutral + blend*(target - neutral);
    end

    %% Selection, event, and publication utilities

    function clearSelection()
        % Clear site and edge selections as one logical operation.
        state.selectedSiteUids = [];
        state.selectedEdgeUids = [];
    end

    function removeUnavailableEdgeSelection()
        % Drop an edge selection after automatic recomputation removes that bond.
        [firstIndex, secondIndex] = selectedEdgeIndices();
        if isempty(firstIndex)
            state.selectedEdgeUids = [];
            return;
        end
        if abs(state.couplings(firstIndex, secondIndex)) <= 1e-14 ...
                && abs(state.couplings(secondIndex, firstIndex)) <= 1e-14 ...
                && ~state.manual(firstIndex, secondIndex) ...
                && ~state.manual(secondIndex, firstIndex)
            state.selectedEdgeUids = [];
        end
    end

    function index = siteIndexFromUid(uid)
        % UIDs remain stable even though display indices are reordered.
        if isempty(uid) || isempty(state.nodes)
            index = [];
            return;
        end
        index = find([state.nodes.uid] == uid, 1, 'first');
    end

    function [firstIndex, secondIndex] = selectedEdgeIndices()
        % Resolve the selected UID pair into current left-to-right site indices.
        firstIndex = [];
        secondIndex = [];
        if numel(state.selectedEdgeUids) ~= 2
            return;
        end
        indexA = siteIndexFromUid(state.selectedEdgeUids(1));
        indexB = siteIndexFromUid(state.selectedEdgeUids(2));
        if isempty(indexA) || isempty(indexB) || indexA == indexB
            return;
        end
        firstIndex = min(indexA, indexB);
        secondIndex = max(indexA, indexB);
        state.selectedEdgeUids = [state.nodes(firstIndex).uid, state.nodes(secondIndex).uid];
    end

    function tf = selectedEdgeMatches(uidA, uidB)
        % Compare an unordered UID pair with the selected edge.
        tf = numel(state.selectedEdgeUids) == 2 ...
            && all(ismember([uidA uidB], state.selectedEdgeUids));
    end

    function position = clickPosition(event)
        % Prefer event intersection data and fall back to the axes current point.
        position = ui.axes.CurrentPoint(1, 1:2);
        if ~isempty(event) && isprop(event, 'IntersectionPoint')
            candidate = event.IntersectionPoint;
            if numel(candidate) >= 2 && all(isfinite(candidate(1:2)))
                position = candidate(1:2);
            end
        end
    end

    function button = mouseButton(event)
        % Normalize sublattice-editor mouse events to MATLAB button codes.
        button = 1;
        if ~isempty(event) && isprop(event, 'Button') && ~isempty(event.Button)
            button = event.Button;
        elseif strcmp(ui.figure.SelectionType, 'alt')
            button = 3;
        end
    end

    %% Hamiltonian assembly and shutdown

    function acceptConfiguration(~, ~)
        % Validate and publish the sublattice model before resuming the caller.
        ui.acceptButton.Enable = 'off';
        drawnow;
        try
            H = calculateHamiltonian();
            configuration = sublatticeConfigurationFromState();
        catch exception
            ui.acceptButton.Enable = 'on';
            uialert(ui.figure, exception.message, 'Unable to build sublattice');
            return;
        end
        if isvalid(ui.figure) && strcmp(ui.figure.WaitStatus, 'waiting')
            uiresume(ui.figure);
        end
    end

    function output = calculateHamiltonian()
        % Transpose (source,destination) storage, then install effective onsite terms.
        count = numel(state.nodes);
        output = effectiveCouplings().';
        if count > 0
            output(1:count + 1:end) = effectiveDetunings();
        end
        if any(~isfinite(output), 'all')
            error('lattice_config:NonfiniteHamiltonian', ...
                'The Hamiltonian contains a nonfinite detuning or coupling.');
        end
    end

    function cancelConfiguration(~, ~)
        % Defer an early close callback until construction can exit safely.
        if ~state.initialized
            state.cancelRequested = true;
            return;
        end
        if isvalid(ui.figure) && strcmp(ui.figure.WaitStatus, 'waiting')
            uiresume(ui.figure);
        elseif isvalid(ui.figure)
            delete(ui.figure);
        end
    end
end
