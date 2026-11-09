function [U,P] = kerr_conservation(plotting,z,c,sysk)
    %% System hierarchy
    narginchk(4,4)

    if ~(isnumeric(c) || (isstruct(c) && isscalar(c) && isfield(c,'full')))
        error('c must be a numeric array or a scalar hierarchy containing c.full.')
    end
    if ~isstruct(sysk) || ~isscalar(sysk)
        error('sysk must be a scalar system structure or hierarchy.')
    end

    [c,c_hierarchical] = system_hierarchy(c);
    [sysk,sysk_hierarchical] = system_hierarchy(sysk);

    if c_hierarchical ~= sysk_hierarchical
        error('c and sysk must use the same leaf or hierarchical convention.')
    end
    if ~isstruct(sysk.full) || ~isscalar(sysk.full)
        error('sysk.full must be a scalar system structure.')
    end

    if sysk_hierarchical
        if ~isfield(sysk,'B') || ~iscell(sysk.B) || isempty(sysk.B) ...
                || ~isstruct(sysk.sub) || ~isscalar(sysk.sub)
            error('A hierarchical sysk must contain packed sub and nonempty B.')
        end
        N_sub = numel(sysk.B);
    else
        N_sub = 0;
    end

    %% Dimensions
    [c.full,N_steps,N_init] = conservation_amplitudes(c.full,z,'c.full');
    M = size(c.full,1);
    z = z(:);

    if numel(sysk.full.epsilon) ~= M
        error('The full-system eigensystem and c must contain the same number of modes.')
    end

    %% Full system
    [U.full,P.full] = kerr_conservation_system(c.full,sysk.full);

    %% Subsystems
    U.sub = repmat(U.full,0,1);
    P.sub = zeros(N_steps,N_init,N_sub);

    if N_sub > 0
        all_indices = [sysk.B{:}];
        if ~isequal(all_indices,1:M)
            error('sysk.B must partition the packed subsystem modes in order.')
        end

        [c.sub,N_sub_steps,N_sub_init] = ...
            conservation_amplitudes(c.sub,z,'c.sub');
        if size(c.sub,1) ~= M
            error('c.sub must contain one packed modal amplitude per full-system mode.')
        end
        if N_sub_steps ~= N_steps || N_sub_init ~= N_init
            error('c.sub must have the same step and member dimensions as c.full.')
        end

        for ss = 1:N_sub
            b = sysk.B{ss};
            subsystem = kerr_conservation_subsystem(sysk.sub,b,M);

            [U.sub(ss),P.sub(:,:,ss)] = ...
                kerr_conservation_system(c.sub(b,:,:),subsystem);
        end
    end

    %% Plot conserved quantities
    if plotting
        kerr_conservation_plots(z,sysk,U,P)
    end

    %% Output convention
    if ~sysk_hierarchical
        U = U.full;
        P = P.full;
    end
end

function [c,N_steps,N_init] = conservation_amplitudes(c,z,name)
    %% Modal-amplitude dimensions
    if ~isnumeric(c) || isempty(c) || ndims(c) > 3
        error('%s must be a nonempty numeric array with at most three dimensions.',name)
    end
    if ~isnumeric(z) || isempty(z) || ~isvector(z)
        error('z must be a nonempty numeric vector.')
    end

    M = size(c,1);
    if ismatrix(c) && isscalar(z)
        N_steps = 1;
        N_init = size(c,2);
    elseif ismatrix(c)
        N_steps = size(c,2);
        N_init = 1;
    else
        N_steps = size(c,2);
        N_init = size(c,3);
    end

    if numel(z) ~= N_steps
        error('z must contain one entry per propagation step in %s.',name)
    end
    c = reshape(c,M,N_steps,N_init);
end

function subsystem = kerr_conservation_subsystem(system,b,M)
    %% Local block from the packed subsystem basis
    I_W = reshape(1:M^2,M,M);
    pair_indices = I_W(b,b);
    pair_indices = pair_indices(:);

    subsystem.epsilon = system.epsilon(b);
    subsystem.G = system.G(b,b);
    subsystem.W = system.W(pair_indices,pair_indices);
    subsystem.K = system.K(b,b);
end

function [U,P] = kerr_conservation_system(c,system)
    %% Dimensions
    M = size(c,1);
    N_steps = size(c,2);
    N_init = size(c,3);

    epsilon = system.epsilon;
    G = system.G;
    W = system.W;
    K = system.K;

    if ~ismatrix(W) || ~isequal(size(W),[M^2 M^2])
        error('W must be an M^2-by-M^2 matrix with pair index m + (n - 1)M.')
    end

    %% Power
    J = abs(c).^2;
    P = zeros(N_steps,N_init);
    for kk = 1:N_init
        cn = c(:,:,kk);
        P(:,kk) = real(sum(conj(cn).*(G*cn),1)).';
    end

    %% Energies
    U = struct;
    U.lin = struct;
    U.non = struct;
    U.tot = struct;
    U.UPA = struct;
    U.FWM = struct;

    %% Exact linear energy
    U.lin.exact = zeros(N_steps,N_init);
    for kk = 1:N_init
        cn = c(:,:,kk);
        U.lin.exact(:,kk) = sum(conj(cn).*(G*(epsilon(:).*cn)),1).';
    end

    %% Exact nonlinear energy
    U.non.exact = zeros(N_steps,N_init);

    N_block = 256;
    for kk = 1:N_init
        for nn = 1:N_block:N_steps
            I_block = nn:min(nn + N_block - 1,N_steps);
            cn = c(:,I_block,kk);

            Q = reshape(conj(cn),M,1,[]).*reshape(cn,1,M,[]);
            Q = reshape(Q,M^2,[]);

            U.non.exact(I_block,kk) = real(sum(Q.*(W*Q),1)).'/2;
        end
    end

    %% Exact total energy
    U.tot.exact = U.lin.exact + U.non.exact;

    %% Uniform phase approximation
    U.UPA.exact = zeros(N_steps,N_init);
    for mm = 1:M
        for nn = 1:M
            J_mn = reshape(J(mm,:,:).*J(nn,:,:),N_steps,N_init);
            U.UPA.exact = U.UPA.exact + K(mm,nn)*J_mn;
        end
    end
    U.UPA.exact = real(U.UPA.exact);

    %% Four-wave mixing
    U.FWM.exact = U.non.exact - U.UPA.exact;

    %% Approximate energies
    U.lin.approx = {U.lin.exact + U.FWM.exact - U.FWM.exact(1,:)};
    U.tot.approx = {U.lin.exact + U.UPA.exact};
end

function kerr_conservation_plots(z,sysk,U,P)
    %% Figure
    viewer = uifigure( ...
        'Name','Kerr conserved quantities', ...
        'Tag','KerrConservationFigure', ...
        'Color',[0.97 0.98 0.99], ...
        'Position',[100 100 1100 700], ...
        'Visible','off');

    viewer_grid = uigridlayout(viewer,[2 1]);
    viewer_grid.RowHeight = {54,'1x'};
    viewer_grid.ColumnWidth = {'1x'};
    viewer_grid.Padding = [8 8 8 8];
    viewer_grid.RowSpacing = 7;

    %% Toolbar
    toolbar = uipanel(viewer_grid, ...
        'BackgroundColor',[1 1 1], ...
        'BorderColor',[0.72 0.72 0.72]);
    toolbar.Layout.Row = 1;
    toolbar.Layout.Column = 1;

    toolbar_grid = uigridlayout(toolbar,[1 4]);
    toolbar_grid.ColumnWidth = {80,180,100,'1x'};
    toolbar_grid.RowHeight = {'1x'};
    toolbar_grid.Padding = [12 7 12 7];
    toolbar_grid.ColumnSpacing = 8;

    quantity_label = uilabel(toolbar_grid, ...
        'Text','Quantity', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center');
    quantity_label.Layout.Column = 1;

    quantity_dropdown = uidropdown(toolbar_grid, ...
        'Items',{'H_lin','H_non','H_tot','H_UPA','H_FWM','P'}, ...
        'ItemsData',{'lin','non','tot','UPA','FWM','power'}, ...
        'Value','lin', ...
        'FontSize',16, ...
        'Tag','KerrQuantityDropdown');
    quantity_dropdown.Layout.Column = 2;

    exact_checkbox = uicheckbox(toolbar_grid, ...
        'Text','Exact', ...
        'Value',true, ...
        'FontSize',16, ...
        'Tag','KerrExactCheckbox');
    exact_checkbox.Layout.Column = 3;

    approximation_grid = uigridlayout(toolbar_grid,[1 1]);
    approximation_grid.Layout.Column = 4;
    approximation_grid.RowHeight = {'1x'};
    approximation_grid.ColumnWidth = {'1x'};
    approximation_grid.Padding = [0 0 0 0];
    approximation_grid.ColumnSpacing = 14;
    approximation_checkboxes = cell(1,0);

    %% Plot panel
    plot_panel = uipanel(viewer_grid, ...
        'BackgroundColor',[1 1 1], ...
        'BorderType','none', ...
        'Tag','KerrConservationPlotPanel');
    plot_panel.Layout.Row = 2;
    plot_panel.Layout.Column = 1;

    %% Callbacks
    quantity_dropdown.ValueChangedFcn = @quantity_changed;
    exact_checkbox.ValueChangedFcn = @selection_changed;

    rebuild_approximation_checkboxes()
    movegui(viewer,'center')
    viewer.Visible = 'on';
    drawnow
    redraw_viewer()

    function selection_changed(~,~)
        redraw_viewer()
    end

    function quantity_changed(~,~)
        rebuild_approximation_checkboxes()
        redraw_viewer()
    end

    function rebuild_approximation_checkboxes()
        delete(approximation_grid.Children)

        quantity = string(quantity_dropdown.Value);
        if quantity ~= "power" && isfield(U.full.(quantity),'approx')
            approximation_names = conservation_approximation_names( ...
                quantity,numel(U.full.(quantity).approx));
        else
            approximation_names = cell(1,0);
        end

        N_approx = numel(approximation_names);
        approximation_checkboxes = cell(1,N_approx);
        if N_approx == 0
            approximation_grid.ColumnWidth = {'1x'};
        else
            approximation_grid.ColumnWidth = repmat({'fit'},1,N_approx);
        end

        for aa = 1:N_approx
            approximation_checkboxes{aa} = uicheckbox(approximation_grid, ...
                'Text',approximation_names{aa}, ...
                'Value',false, ...
                'FontSize',16, ...
                'Tag',sprintf('KerrApproxCheckbox%d',aa));
            approximation_checkboxes{aa}.Layout.Column = aa;
            approximation_checkboxes{aa}.ValueChangedFcn = @selection_changed;
        end
    end

    function redraw_viewer()
        if ~isvalid(viewer) || ~isvalid(plot_panel)
            return
        end

        quantity = string(quantity_dropdown.Value);
        has_approximation = ~isempty(approximation_checkboxes);

        if quantity ~= "power" && ~has_approximation && ~exact_checkbox.Value
            exact_checkbox.Value = true;
        end

        quantity_dropdown.Enable = 'off';
        exact_checkbox.Enable = 'off';
        for aa = 1:numel(approximation_checkboxes)
            approximation_checkboxes{aa}.Enable = 'off';
        end
        viewer.Pointer = 'watch';
        drawnow

        cleanup = onCleanup(@() restore_controls( ...
            viewer,quantity_dropdown,exact_checkbox,approximation_checkboxes));

        show_approx = false(1,numel(approximation_checkboxes));
        for aa = 1:numel(approximation_checkboxes)
            show_approx(aa) = approximation_checkboxes{aa}.Value;
        end

        delete(plot_panel.Children)
        plot_conservation_quantity(plot_panel,z,sysk,U,P,quantity, ...
            exact_checkbox.Value,show_approx)

        clear cleanup
    end
end

function restore_controls(viewer,quantity_dropdown,exact_checkbox,approximation_checkboxes)
    if isvalid(viewer)
        viewer.Pointer = 'arrow';
    end
    if isvalid(quantity_dropdown)
        quantity_dropdown.Enable = 'on';
    end
    if ~isvalid(exact_checkbox)
        return
    end

    quantity = string(quantity_dropdown.Value);
    if quantity == "power"
        exact_checkbox.Enable = 'off';
    else
        exact_checkbox.Enable = 'on';
    end

    for aa = 1:numel(approximation_checkboxes)
        if isvalid(approximation_checkboxes{aa})
            approximation_checkboxes{aa}.Enable = 'on';
        end
    end
end

function plot_conservation_quantity(parent,z,sysk,U,P,quantity,show_exact,show_approx)
    %% Axes
    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    hold(plot_axes,'on')

    %% Systems
    N_sub = numel(U.sub);
    if N_sub == 1
        N_sub = 0;
    end
    N_systems = N_sub + 1;
    [system_colors_start,system_colors_end] = conservation_system_colors(N_sub);
    member_opacity = 0.2;
    energy_is_complex = false;
    member_marker = 'none';
    mean_marker = 'none';
    if isscalar(z)
        member_marker = '.';
        mean_marker = 'o';
    end

    if quantity == "power"
        N_signals = N_systems;
    else
        N_signals_per_system = double(show_exact);
        if any(show_approx) && isfield(U.full.(quantity),'approx')
            N_signals_per_system = N_signals_per_system ...
                + nnz(show_approx);
        end
        N_signals = N_systems*N_signals_per_system;
    end

    signal_template = struct('members',[],'color',[0 0 0],'name','');
    signals = repmat(signal_template,N_signals,1);
    signal_index = 0;

    for ss = 1:N_systems
        if ss == 1
            system_name = 'Full system';
            U_system = U.full;
            P_system = P.full;
        else
            system_name = sprintf('Subsystem %d',ss - 1);
            U_system = U.sub(ss - 1);
            P_system = P.sub(:,:,ss - 1);
        end

        line_color = system_colors_start(ss,:);

        if quantity == "power"
            signal_index = signal_index + 1;
            signals(signal_index).members = P_system;
            signals(signal_index).color = line_color;
            signals(signal_index).name = system_name;
            continue
        end

        component = U_system.(quantity);
        if show_exact
            signal_index = signal_index + 1;
            signals(signal_index).members = component.exact;
            signals(signal_index).color = line_color;
            signals(signal_index).name = sprintf('%s: Exact',system_name);
            energy_is_complex = energy_is_complex || has_imaginary_part(component.exact);
        end

        if any(show_approx) && isfield(component,'approx')
            approximation_names = conservation_approximation_names(quantity,numel(component.approx));
            for aa = 1:numel(component.approx)
                if aa > numel(show_approx) || ~show_approx(aa)
                    continue
                end

                approximation = component.approx{aa};
                fraction = aa/numel(component.approx);
                approximation_color = (1 - fraction)*system_colors_start(ss,:) ...
                    + fraction*system_colors_end(ss,:);

                signal_index = signal_index + 1;
                signals(signal_index).members = approximation;
                signals(signal_index).color = approximation_color;
                signals(signal_index).name = ...
                    sprintf('%s: %s',system_name,approximation_names{aa});
                energy_is_complex = energy_is_complex || has_imaginary_part(approximation);
            end
        end
    end

    %% Member signals
    y_min = inf;
    y_max = -inf;
    for qq = 1:N_signals
        signal = real(signals(qq).members);
        member_color = member_opacity*signals(qq).color + (1 - member_opacity);
        z_members = repmat(z,1,size(signal,2));

        plot(plot_axes,z_members,signal,':','Color',member_color, ...
            'Marker',member_marker,'HandleVisibility','off')

        finite_signal = signal(isfinite(signal));
        if ~isempty(finite_signal)
            y_min = min(y_min,min(finite_signal));
            y_max = max(y_max,max(finite_signal));
        end
    end

    %% Ensemble averages
    for qq = 1:N_signals
        plot(plot_axes,z,mean(real(signals(qq).members),2),'-', ...
            'Color',signals(qq).color,'LineWidth',2,'Marker',mean_marker, ...
            'DisplayName',signals(qq).name)
    end

    %% Labels
    hold(plot_axes,'off')
    grid(plot_axes,'on')
    box(plot_axes,'on')
    xlabel(plot_axes,'Propagation coordinate')

    quantity_label = conservation_quantity_label(quantity,energy_is_complex);
    title(plot_axes,quantity_label,'Interpreter','latex')
    ylabel(plot_axes,quantity_label,'Interpreter','latex')

    if ~isempty(plot_axes.Children)
        legend(plot_axes,'Location','best')
    end

    %% Linear-energy limits
    if quantity == "lin"
        epsilon = sysk.full.epsilon;
        G = sysk.full.G;
        P0 = mean(P.full(1,:));
        epsilon_scale = max(1,max(abs(epsilon),[],'all'));
        G_scale = max(1,norm(G,'fro'));
        G_is_identity = norm(G - eye(size(G)),'fro') <= 1e-10*G_scale;
        if all(abs(imag(epsilon)) <= 1e-12*epsilon_scale,'all') && G_is_identity
            H_min = min(real(epsilon),[],'all')*P0;
            H_max = max(real(epsilon),[],'all')*P0;
            hold(plot_axes,'on')
            yline(plot_axes,H_min,'-', ...
                '$\min(H_{\mathrm{lin}})=\varepsilon_1P$', ...
                'Interpreter','latex','LabelVerticalAlignment','bottom', ...
                'HandleVisibility','off')
            yline(plot_axes,H_max,'-', ...
                '$\max(H_{\mathrm{lin}})=\varepsilon_MP$', ...
                'Interpreter','latex','LabelVerticalAlignment','top', ...
                'HandleVisibility','off')
            hold(plot_axes,'off')

            y_min = min(y_min,H_min);
            y_max = max(y_max,H_max);
        end
    end

    %% Vertical limits
    if isfinite(y_min) && isfinite(y_max)
        y_span = y_max - y_min;
        y_scale = max(abs([y_min y_max]));
        y_padding = max(0.08*y_span,0.01*y_scale);
        if y_padding == 0
            y_padding = 1;
        end
        y_lower = y_min - y_padding;
        if quantity == "power"
            y_lower = max(0,y_lower);
        end
        ylim(plot_axes,[y_lower y_max + y_padding])
    end
end

function names = conservation_approximation_names(quantity,N_approx)
    switch quantity
        case "lin"
            names = {'FWM-shifted'};
        case "tot"
            names = {'Linear + UPA'};
        otherwise
            names = cell(1,0);
    end

    for aa = numel(names) + 1:N_approx
        names{aa} = sprintf('Approximation %d',aa);
    end
end

function [start_colors,end_colors] = conservation_system_colors(N_sub)
    %% Full system
    start_colors = zeros(N_sub + 1,3);
    end_colors = zeros(N_sub + 1,3);
    start_colors(1,:) = [0 0 0];
    end_colors(1,:) = [0.65 0.65 0.65];

    %% Subsystems
    start_hues = [2/3 0 1/3 5/6 1/2 1/6];
    for ss = 1:N_sub
        hue_index = 1 + mod(ss - 1,numel(start_hues));
        color_cycle = floor((ss - 1)/numel(start_hues));
        saturation = max(0.55,1 - 0.15*color_cycle);
        start_value = max(0.55,0.85 - 0.08*color_cycle);
        end_value = max(0.65,0.95 - 0.08*color_cycle);

        start_hue = start_hues(hue_index);
        end_hue = mod(start_hue - 1/6,1);
        start_colors(ss + 1,:) = hsv2rgb([start_hue saturation start_value]);
        end_colors(ss + 1,:) = hsv2rgb([end_hue saturation end_value]);
    end
end

function label = conservation_quantity_label(quantity,energy_is_complex)
    switch quantity
        case "lin"
            label = 'H_{\mathrm{lin}}';
        case "non"
            label = 'H_{\mathrm{non}}';
        case "tot"
            label = 'H_{\mathrm{tot}}';
        case "UPA"
            label = 'H_{\mathrm{UPA}}';
        case "FWM"
            label = 'H_{\mathrm{FWM}}';
        case "power"
            label = 'P';
    end

    if energy_is_complex
        label = append('$\mathrm{Re}\,',label,'$');
    else
        label = append('$',label,'$');
    end
end

function value = has_imaginary_part(signal)
    signal_scale = max(1,max(abs(signal),[],'all'));
    value = any(abs(imag(signal)) > 1e-12*signal_scale,'all');
end
