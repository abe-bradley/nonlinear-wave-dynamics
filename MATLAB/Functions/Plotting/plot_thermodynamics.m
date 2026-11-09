function plot_thermodynamics(statistics,U,P,sysk)
    %% System hierarchy
    narginchk(4,4)

    [statistics,statistics_hierarchical] = system_hierarchy(statistics);
    [U,U_hierarchical] = system_hierarchy(U);
    [P,P_hierarchical] = system_hierarchy(P);
    [sysk,sysk_hierarchical] = system_hierarchy(sysk);

    hierarchy = [statistics_hierarchical,U_hierarchical, ...
        P_hierarchical,sysk_hierarchical];
    if any(hierarchy ~= hierarchy(1))
        error('statistics, U, P, and sysk must use the same system hierarchy.')
    end
    hierarchical = hierarchy(1);

    if ~isscalar(statistics) || ~isscalar(U) || ~isscalar(P) ...
            || ~isscalar(sysk) || ~isstruct(statistics.full) ...
            || ~isscalar(statistics.full) || ~isstruct(U.full) ...
            || ~isscalar(U.full) || ~isnumeric(P.full) ...
            || isempty(P.full) || ~ismatrix(P.full) ...
            || ~isstruct(sysk.full) || ~isscalar(sysk.full)
        error('Inputs must be scalar system containers with valid full fields.')
    end

    %% Systems
    M = thermodynamic_mode_count(statistics.full,'statistics.full');
    if ~isfield(sysk.full,'epsilon') || numel(sysk.full.epsilon) ~= M
        error('sysk.full must contain one energy per full-system mode.')
    end

    B = cell(0,1);
    N_sub = 0;
    if hierarchical
        if ~isfield(sysk,'B') || ~iscell(sysk.B) || isempty(sysk.B) ...
                || ~isstruct(sysk.sub) || ~isscalar(sysk.sub)
            error('A hierarchical sysk must contain packed sub and nonempty B.')
        end

        B = sysk.B(:).';
        N_sub = numel(B);
        if ~isstruct(statistics.sub) || numel(statistics.sub) ~= N_sub ...
                || ~isstruct(U.sub) || numel(U.sub) ~= N_sub ...
                || ~isnumeric(P.sub) || isempty(P.sub) || ndims(P.sub) > 3 ...
                || size(P.sub,1) ~= size(P.full,1) ...
                || size(P.sub,2) ~= size(P.full,2) ...
                || size(P.sub,3) ~= N_sub
            error('The subsystem statistics, U, and P must match sysk.B.')
        end
        if ~isfield(sysk.sub,'epsilon') || numel(sysk.sub.epsilon) ~= M
            error('sysk.sub must contain one packed energy per full-system mode.')
        end

        all_indices = zeros(1,0);
        for sub_index = 1:N_sub
            b = B{sub_index};
            if ~isnumeric(b) || ~isreal(b) || isempty(b) || ~isvector(b) ...
                    || any(~isfinite(b)) || any(b ~= round(b)) ...
                    || any(b < 1) || any(b > M)
                error('sysk.B{%d} contains invalid modal indices.',sub_index)
            end
            B{sub_index} = b(:).';
            all_indices = [all_indices B{sub_index}]; %#ok<AGROW>

            M_sub = thermodynamic_mode_count( ...
                statistics.sub(sub_index), ...
                sprintf('statistics.sub(%d)',sub_index));
            if M_sub ~= numel(B{sub_index})
                error('statistics.sub(%d) does not match sysk.B{%d}.', ...
                    sub_index,sub_index)
            end
        end
        if ~isequal(all_indices,1:M)
            error('sysk.B must partition the packed subsystem modes in order.')
        end
    end

    if N_sub <= 1
        N_sub_plot = 0;
    else
        N_sub_plot = N_sub;
    end

    %% Figure
    viewer = uifigure( ...
        'Name','Thermodynamics', ...
        'Tag','ThermodynamicsFigure', ...
        'Color',[0.97 0.98 0.99], ...
        'Visible','off');
    viewer.WindowState = 'maximized';

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

    toolbar_grid = uigridlayout(toolbar,[1 7]);
    toolbar_grid.ColumnWidth = {70,220,55,250,55,100,'1x'};
    toolbar_grid.RowHeight = {'1x'};
    toolbar_grid.Padding = [12 7 12 7];
    toolbar_grid.ColumnSpacing = 8;

    system_label = uilabel(toolbar_grid, ...
        'Text','System', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center');
    system_label.Layout.Column = 1;

    system_names = [{'Full system'},arrayfun( ...
        @(ss) sprintf('Subsystem %d',ss),1:N_sub_plot,'UniformOutput',false)];
    system_dropdown = uidropdown(toolbar_grid, ...
        'Items',system_names, ...
        'ItemsData',1:numel(system_names), ...
        'Value',1, ...
        'FontSize',16, ...
        'Tag','ThermodynamicsSystemDropdown');
    system_dropdown.Layout.Column = 2;

    plot_label = uilabel(toolbar_grid, ...
        'Text','Plot', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center');
    plot_label.Layout.Column = 3;

    plot_dropdown = uidropdown(toolbar_grid, ...
        'Items',{'Mean occupations','Reciprocal occupations', ...
        'Marginals','Pearson correlations','Two-mode joint'}, ...
        'ItemsData',{'occupations','reciprocal','marginals','pearson','joint'}, ...
        'Value','occupations', ...
        'FontSize',16, ...
        'Tag','ThermodynamicsPlotDropdown');
    plot_dropdown.Layout.Column = 4;

    mode_label = uilabel(toolbar_grid, ...
        'Text','Mode', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center', ...
        'Visible','off', ...
        'Tag','ThermodynamicsModeLabel');
    mode_label.Layout.Column = 5;

    mode_dropdown = uidropdown(toolbar_grid, ...
        'Items',{'1'}, ...
        'Value','1', ...
        'FontSize',16, ...
        'Visible','off', ...
        'Tag','ThermodynamicsModeDropdown');
    mode_dropdown.Layout.Column = 6;

    %% Plot panel
    plot_panel = uipanel(viewer_grid, ...
        'BackgroundColor',[1 1 1], ...
        'BorderType','none', ...
        'Tag','ThermodynamicsPlotPanel');
    plot_panel.Layout.Row = 2;
    plot_panel.Layout.Column = 1;

    %% Callbacks
    system_dropdown.ValueChangedFcn = @system_changed;
    plot_dropdown.ValueChangedFcn = @plot_changed;
    mode_dropdown.ValueChangedFcn = @selection_changed;

    update_mode_selector()
    update_mode_visibility()
    viewer.Visible = 'on';
    drawnow
    redraw_viewer()

    function system_changed(~,~)
        update_mode_selector()
        redraw_viewer()
    end

    function plot_changed(~,~)
        update_mode_visibility()
        redraw_viewer()
    end

    function selection_changed(~,~)
        redraw_viewer()
    end

    function update_mode_selector()
        [statistics_plot,~,~,~] = selected_system();
        N_modes = thermodynamic_mode_count(statistics_plot,'statistics');
        old_mode = str2double(string(mode_dropdown.Value));
        if ~isfinite(old_mode)
            old_mode = 1;
        end

        mode_items = arrayfun(@num2str,1:N_modes,'UniformOutput',false);
        mode_dropdown.Items = mode_items;
        mode_dropdown.Value = mode_items{min(max(round(old_mode),1),N_modes)};
    end

    function update_mode_visibility()
        visible = 'off';
        if string(plot_dropdown.Value) == "joint"
            visible = 'on';
        end
        mode_label.Visible = visible;
        mode_dropdown.Visible = visible;
    end

    function [statistics_plot,U_plot,P_plot,epsilon_plot] = selected_system()
        system_index = system_dropdown.Value;
        if system_index == 1
            statistics_plot = statistics.full;
            U_plot = U.full;
            P_plot = P.full;
            epsilon_plot = sysk.full.epsilon;
        else
            selected_subsystem = system_index - 1;
            statistics_plot = statistics.sub(selected_subsystem);
            U_plot = U.sub(selected_subsystem);
            P_plot = P.sub(:,:,selected_subsystem);
            epsilon_plot = sysk.sub.epsilon(B{selected_subsystem});
        end
    end

    function redraw_viewer()
        if ~isvalid(viewer) || ~isvalid(plot_panel)
            return
        end

        system_dropdown.Enable = 'off';
        plot_dropdown.Enable = 'off';
        mode_dropdown.Enable = 'off';
        viewer.Pointer = 'watch';
        drawnow
        cleanup = onCleanup(@() thermodynamic_restore_controls( ...
            viewer,system_dropdown,plot_dropdown,mode_dropdown));

        delete(plot_panel.Children)
        try
            [statistics_plot,U_plot,P_plot,epsilon_plot] = selected_system();
            plot_type = string(plot_dropdown.Value);

            if plot_type == "joint"
                mode_index = str2double(string(mode_dropdown.Value));
                plot_two_mode_joint(plot_panel,statistics_plot,P_plot,mode_index)
            else
                [U_state,P_state,epsilon_plot] = thermodynamic_state( ...
                    statistics_plot,U_plot,P_plot,epsilon_plot);

                switch plot_type
                    case "occupations"
                        plot_mean_occupations(plot_panel,statistics_plot, ...
                            U_state,P_state,epsilon_plot)
                    case "reciprocal"
                        plot_reciprocal_occupations(plot_panel,statistics_plot, ...
                            U_state,P_state,epsilon_plot)
                    case "marginals"
                        plot_marginals(plot_panel,statistics_plot, ...
                            U_state,P_state,epsilon_plot)
                    case "pearson"
                        plot_pearson(plot_panel,statistics_plot, ...
                            U_state,P_state,epsilon_plot)
                end
            end
        catch exception
            delete(plot_panel.Children)
            thermodynamic_message(plot_panel,exception.message)
        end

        drawnow
        clear cleanup
    end
end

function M = thermodynamic_mode_count(statistics,name)
    if ~isstruct(statistics) || ~isscalar(statistics) ...
            || ~isfield(statistics,'N_modes') ...
            || ~isscalar(statistics.N_modes) || statistics.N_modes < 1 ...
            || statistics.N_modes ~= round(statistics.N_modes)
        error('%s must contain a positive integer N_modes.',name)
    end
    M = statistics.N_modes;

    if ~isfield(statistics,'first_moment') ...
            || numel(statistics.first_moment) ~= M
        error('%s.first_moment must contain one value per mode.',name)
    end
end

function [U_state,P_state,epsilon] = ...
        thermodynamic_state(statistics,U,P,epsilon)
    %% Late-window state
    M = thermodynamic_mode_count(statistics,'statistics');
    if M < 3
        error('Thermodynamic predictions require at least three modes.')
    end
    if ~isfield(statistics,'N_steps') || ~isfield(statistics,'N_init')
        error('statistics must contain N_steps and N_init.')
    end
    if ~isstruct(U) || ~isscalar(U) || ~isfield(U,'lin') ...
            || ~isfield(U.lin,'exact')
        error('U must contain the exact linear energy in U.lin.exact.')
    end
    if ~isnumeric(P) || isempty(P) || ~ismatrix(P)
        error('P must contain the step-by-member power signal.')
    end

    U_lin = U.lin.exact;
    N_steps = statistics.N_steps;
    N_init = statistics.N_init;
    if ~isnumeric(U_lin) || isempty(U_lin) || ~ismatrix(U_lin) ...
            || size(U_lin,1) < N_steps || size(P,1) < N_steps ...
            || size(U_lin,2) ~= N_init || size(P,2) ~= N_init
        error('U and P must cover the members and late steps used in statistics.')
    end

    I_end_U = size(U_lin,1) - N_steps + 1:size(U_lin,1);
    I_end_P = size(P,1) - N_steps + 1:size(P,1);
    U_state = mean(U_lin(I_end_U,:),'all');
    P_state = mean(P(I_end_P,:),'all');
    epsilon = epsilon(:);

    energy_scale = max(1,max(abs([epsilon;U_state])));
    if max(abs(imag([epsilon;U_state]))) > 1e-10*energy_scale
        error('Thermodynamic predictions require real modal energies and linear energy.')
    end
    epsilon = real(epsilon);
    U_state = real(U_state);

    if numel(epsilon) ~= M || any(~isfinite(epsilon)) ...
            || ~isfinite(U_state)
        error('The modal energies and mean linear energy must be finite.')
    end
    if ~isreal(P_state) || ~isfinite(P_state) || P_state <= 0
        error('The mean late-window power must be finite, real, and positive.')
    end
end

function thermodynamic_restore_controls( ...
        viewer,system_dropdown,plot_dropdown,mode_dropdown)
    if isvalid(viewer)
        viewer.Pointer = 'arrow';
    end
    if isvalid(system_dropdown)
        system_dropdown.Enable = 'on';
    end
    if isvalid(plot_dropdown)
        plot_dropdown.Enable = 'on';
    end
    if isvalid(mode_dropdown)
        mode_dropdown.Enable = 'on';
    end
end

function thermodynamic_message(parent,message)
    message_grid = uigridlayout(parent,[1 1]);
    message_grid.Padding = [0 0 0 0];
    message_label = uilabel(message_grid, ...
        'Text',message, ...
        'FontSize',20, ...
        'FontColor',[0.45 0.48 0.52], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','center', ...
        'WordWrap','on');
    message_label.Layout.Row = 1;
    message_label.Layout.Column = 1;
end

function plot_mean_occupations(parent,statistics,U,P,epsilon)
    %% Thermodynamic predictions
    [~,~,J_RJ] = calculate_RJ(U,P,epsilon);
    [~,~,J_MC] = calculate_MC(U,P,epsilon);

    M = statistics.N_modes;
    modes = (1:M).';
    J_sim = statistics.first_moment(:);

    %% Plot
    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    hold(plot_axes,'on')

    mango_stops = [
        0.96 0.25 0.02;
        1.00 0.58 0.00;
        1.00 0.86 0.08
    ];
    mango_colors = interp1([0 0.5 1],mango_stops,linspace(0,0.8,64));
    for mm = 1:M
        gradient_bar(plot_axes,mm,J_sim(mm)/P,0.6,mango_colors);
    end

    simulation_plot = patch(plot_axes,nan,nan,[1.00 0.50 0.00], ...
        'FaceAlpha',0.92,'EdgeColor','none');

    marker_offset = 0.09;
    RJ_modes = modes - marker_offset;
    MC_modes = modes + marker_offset;
    interpolation_modes = linspace(1,M,max(200,20*M));
    RJ_interpolation = interp1(RJ_modes,J_RJ/P, ...
        interpolation_modes - marker_offset,'pchip');
    MC_interpolation = interp1(MC_modes,J_MC/P, ...
        interpolation_modes + marker_offset,'pchip');

    plot(plot_axes,interpolation_modes - marker_offset,RJ_interpolation, ...
        ':','Color',[0.32 0.32 0.32],'LineWidth',1.25, ...
        'HandleVisibility','off')
    plot(plot_axes,RJ_modes,J_RJ/P,'d','LineStyle','none', ...
        'Color',[0.32 0.32 0.32],'MarkerFaceColor','w', ...
        'LineWidth',1.25,'MarkerSize',7.2,'HandleVisibility','off')
    plot(plot_axes,interpolation_modes + marker_offset,MC_interpolation, ...
        ':','Color',[0 0 0],'LineWidth',1.25,'HandleVisibility','off')
    plot(plot_axes,MC_modes,J_MC/P,'o','LineStyle','none', ...
        'Color',[0 0 0],'MarkerFaceColor',[0 0 0], ...
        'LineWidth',1.25,'MarkerSize',7.8,'HandleVisibility','off')

    RJ_plot = plot(plot_axes,nan,nan,':d', ...
        'Color',[0.32 0.32 0.32],'MarkerFaceColor','w', ...
        'LineWidth',1.25,'MarkerSize',7.2);
    MC_plot = plot(plot_axes,nan,nan,':o', ...
        'Color',[0 0 0],'MarkerFaceColor',[0 0 0], ...
        'LineWidth',1.25,'MarkerSize',7.8);

    hold(plot_axes,'off')
    xlabel(plot_axes,'Mode index, m')
    ylabel(plot_axes,'\langle J_m\rangle/P')
    title(plot_axes,'Equilibrium mean occupations')
    legend(plot_axes,[simulation_plot,RJ_plot,MC_plot], ...
        {'Simulated','Rayleigh-Jeans','Microcanonical'}, ...
        'Location','northeast','Box','off')
    xlim(plot_axes,[0.35 M + 0.65])
    xticks(plot_axes,modes)

    occupation_max = max([J_sim;J_RJ(:);J_MC(:)])/P;
    if ~isfinite(occupation_max) || occupation_max <= 0
        occupation_max = 1;
    end
    ylim(plot_axes,[0 1.1*occupation_max])
    grid(plot_axes,'on')
    box(plot_axes,'off')
    plot_axes.FontName = 'Times New Roman';
    plot_axes.FontSize = 14;
    plot_axes.LineWidth = 1.1;
    plot_axes.GridAlpha = 0.14;
    plot_axes.YMinorGrid = 'on';
    plot_axes.TickDir = 'out';
    plot_axes.Layer = 'top';
end

function plot_reciprocal_occupations(parent,statistics,U,P,epsilon)
    %% Thermodynamic predictions
    [~,~,J_RJ] = calculate_RJ(U,P,epsilon);
    [~,~,J_MC] = calculate_MC(U,P,epsilon);

    J_sim = statistics.first_moment(:);
    reciprocal_sim = reciprocal_occupation(J_sim/P);
    reciprocal_RJ = reciprocal_occupation(J_RJ/P);
    reciprocal_MC = reciprocal_occupation(J_MC/P);
    [epsilon_RJ,reciprocal_RJ_plot] = ...
        connect_energy_bands(epsilon,reciprocal_RJ);
    [epsilon_MC,reciprocal_MC_plot] = ...
        connect_energy_bands(epsilon,reciprocal_MC);

    %% Plot
    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    hold(plot_axes,'on')

    simulated_plot = plot(plot_axes,epsilon,reciprocal_sim,'o', ...
        'LineStyle','none','Color',[1.00 0.50 0.00], ...
        'MarkerFaceColor',[1.00 0.50 0.00],'MarkerSize',5);
    RJ_plot = plot(plot_axes,epsilon_RJ,reciprocal_RJ_plot,':', ...
        'Color',[0.20 0.20 0.20],'LineWidth',2);
    MC_plot = plot(plot_axes,epsilon_MC,reciprocal_MC_plot,'-', ...
        'Color',[0.50 0.10 0.40],'LineWidth',2);

    plot(plot_axes,epsilon,reciprocal_RJ,'d','LineStyle','none', ...
        'Color',[0.20 0.20 0.20],'MarkerFaceColor','w', ...
        'MarkerSize',5,'HandleVisibility','off')
    plot(plot_axes,epsilon,reciprocal_MC,'o','LineStyle','none', ...
        'Color',[0.50 0.10 0.40],'MarkerFaceColor',[0.50 0.10 0.40], ...
        'MarkerSize',4,'HandleVisibility','off')

    hold(plot_axes,'off')
    xlabel(plot_axes,'Energy, \epsilon_m')
    ylabel(plot_axes,'P/\langle J_m\rangle')
    title(plot_axes,'Reciprocal occupations')
    legend(plot_axes,[simulated_plot,RJ_plot,MC_plot], ...
        {'Simulated','Rayleigh-Jeans','Microcanonical'}, ...
        'Location','eastoutside','Box','off')
    xlim(plot_axes,tight_plot_limits(epsilon))

    reciprocal_values = [reciprocal_sim;reciprocal_RJ;reciprocal_MC];
    reciprocal_values = reciprocal_values(isfinite(reciprocal_values));
    if ~isempty(reciprocal_values)
        ylim(plot_axes,tight_plot_limits(reciprocal_values))
    end
    grid(plot_axes,'on')
    box(plot_axes,'off')
    plot_axes.TickDir = 'out';
    plot_axes.Layer = 'top';
end

function plot_marginals(parent,statistics,U,P,epsilon)
    %% Plot
    M = statistics.N_modes;
    plot_layout = tiledlayout(parent,'flow', ...
        'TileSpacing','compact','Padding','compact');
    marginal_x = linspace(0,P,200);
    grape_bottom = [0.66 0.75 1.00];
    grape_top = [0.13 0.20 0.52];

    for mm = 1:M
        plot_axes = nexttile(plot_layout);
        hold(plot_axes,'on')

        marginal_sim_edges = statistics.marginal.edges(mm,:)/P;
        marginal_sim_pdf = P*statistics.marginal.pdf(mm,:);
        [marginal_MC_x,marginal_MC_pdf] = ...
            marginals(mm,marginal_x,U,P,epsilon);

        simulation_plot = gradient_histogram(plot_axes, ...
            marginal_sim_edges,marginal_sim_pdf,grape_bottom,grape_top);
        theory_plot = plot(plot_axes,marginal_MC_x/P, ...
            P*marginal_MC_pdf,'Color',[1.00 0.40 0.60], ...
            'LineWidth',2.8);

        hold(plot_axes,'off')
        title(plot_axes,append("m = ",string(mm)))
        xlim(plot_axes,[0 1])

        marginal_max = max([marginal_sim_pdf,marginal_MC_pdf*P],[],'all');
        if ~isfinite(marginal_max) || marginal_max <= 0
            marginal_max = 1;
        end
        ylim(plot_axes,[0 1.1*marginal_max])
        grid(plot_axes,'on')
        box(plot_axes,'on')

        if mm == 1
            legend(plot_axes,[simulation_plot,theory_plot], ...
                {'Simulated','Microcanonical'},'Location','best','Box','off')
        end
    end

    title(plot_layout,'Marginal distributions')
    xlabel(plot_layout,'J_m/P')
    ylabel(plot_layout,'p_m(J_m/P)')
end

function plot_pearson(parent,statistics,U,P,epsilon)
    %% Correlations
    M = statistics.N_modes;
    modes = (1:M).';
    pearson_sim = statistics.pearson_correlations;
    pearson_MC = pearson_correlations(U,P,epsilon);

    pearson_sim(~isfinite(pearson_sim)) = 0;
    pearson_sim = max(min(real(pearson_sim),1),-1);
    pearson_MC(~isfinite(pearson_MC)) = 0;
    pearson_MC = max(min(real(pearson_MC),1),-1);
    pearson_sim = triu(pearson_sim,1);
    pearson_MC = triu(pearson_MC,1);

    correlation_colors = flipud(ot_colormap('redshift',256));
    color_values = linspace(-1,1,size(correlation_colors,1)).';
    stretched_values = sign(color_values).*abs(color_values).^(1/3);
    correlation_colors = interp1( ...
        color_values,correlation_colors,stretched_values,'linear');

    %% Plot
    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    pearson_bars = bar3(plot_axes,abs(pearson_sim));
    hold(plot_axes,'on')

    for nn = 1:numel(pearson_bars)
        correlation_sign = repmat( ...
            repelem(sign(pearson_sim(:,nn)),6),1, ...
            size(pearson_bars(nn).ZData,2));
        pearson_bars(nn).CData = correlation_sign.*pearson_bars(nn).ZData;

        zero_rows = repelem(pearson_sim(:,nn) == 0,6);
        pearson_bars(nn).XData(zero_rows,:) = nan;
        pearson_bars(nn).YData(zero_rows,:) = nan;
        pearson_bars(nn).ZData(zero_rows,:) = nan;
        pearson_bars(nn).CData(zero_rows,:) = nan;
        pearson_bars(nn).FaceColor = 'interp';
        pearson_bars(nn).FaceAlpha = 0.5;
        pearson_bars(nn).EdgeColor = [0 0 0];
        pearson_bars(nn).LineWidth = 0.15;
    end

    [mode_n,mode_m] = meshgrid(modes,modes);
    theory_points = triu(true(M),1) & pearson_MC ~= 0;
    pearson_theory = scatter3(plot_axes,mode_n(theory_points), ...
        mode_m(theory_points),abs(pearson_MC(theory_points)),42, ...
        pearson_MC(theory_points),'filled','MarkerEdgeColor',[0.05 0.05 0.05], ...
        'LineWidth',0.8);

    hold(plot_axes,'off')
    colormap(plot_axes,correlation_colors)
    clim(plot_axes,[-1 1])
    pearson_colorbar = colorbar(plot_axes);
    pearson_colorbar.Label.String = 'Pearson correlation';
    xlabel(plot_axes,'m')
    ylabel(plot_axes,'n')
    zlabel(plot_axes,'|\rho_{mn}|')
    title(plot_axes,'Pearson correlations')
    xlim(plot_axes,[0.5 M + 0.5])
    ylim(plot_axes,[0.5 M + 0.5])
    zlim(plot_axes,[0 1])
    view(plot_axes,[-30 25])
    grid(plot_axes,'on')
    box(plot_axes,'off')
    triangular_prism(plot_axes,[0.5 M + 0.5],[0 1])

    pearson_bars(1).DisplayName = 'Simulated';
    pearson_theory.DisplayName = 'Microcanonical';
    legend(plot_axes,[pearson_bars(1),pearson_theory], ...
        'Location','best','Box','off')
end

function plot_two_mode_joint(parent,statistics,P,mode_index)
    %% Joint distributions
    M = statistics.N_modes;
    if M < 2 || ~isfield(statistics,'joint') ...
            || statistics.joint.N_pairs < 1
        thermodynamic_message(parent,'No two-mode joint distributions are available.')
        return
    end
    if ~isscalar(mode_index) || mode_index < 1 || mode_index > M ...
            || mode_index ~= round(mode_index)
        error('The selected mode index is invalid.')
    end

    N_steps = statistics.N_steps;
    N_init = statistics.N_init;
    if ~isnumeric(P) || ~ismatrix(P) || size(P,1) < N_steps ...
            || size(P,2) ~= N_init
        error('P must cover the members and late steps used in statistics.')
    end
    I_end = size(P,1) - N_steps + 1:size(P,1);
    P_state = mean(P(I_end,:),'all');
    if ~isreal(P_state) || ~isfinite(P_state) || P_state <= 0
        error('The mean late-window power must be finite, real, and positive.')
    end

    modes = (1:M).';
    other_modes = modes(modes ~= mode_index);
    plot_layout = tiledlayout(parent,'flow', ...
        'TileSpacing','compact','Padding','compact');
    purple_colors = ot_colormap('Purples').^2;

    for kk = 1:numel(other_modes)
        other_mode = other_modes(kk);
        mode_pair = sort([mode_index other_mode]);
        pair_index = find(all( ...
            statistics.joint.mode_pairs == mode_pair,2),1);
        if isempty(pair_index)
            error('The stored statistics do not contain modes %d and %d.', ...
                mode_index,other_mode)
        end

        joint_x = statistics.joint.centers(mode_index,:)/P_state;
        joint_y = statistics.joint.centers(other_mode,:)/P_state;
        if mode_index == mode_pair(1)
            joint_pdf = P_state^2*statistics.joint.pdf(:,:,pair_index).';
        else
            joint_pdf = P_state^2*statistics.joint.pdf(:,:,pair_index);
        end

        plot_axes = nexttile(plot_layout);
        imagesc(plot_axes,joint_x,joint_y,joint_pdf)
        plot_axes.YDir = 'normal';
        colormap(plot_axes,purple_colors)

        selected_edges = statistics.joint.edges(mode_index,:)/P_state;
        other_edges = statistics.joint.edges(other_mode,:)/P_state;
        xlim(plot_axes,[max(0,selected_edges(1)) selected_edges(end)])
        ylim(plot_axes,[max(0,other_edges(1)) other_edges(end)])
        title(plot_axes,append("n = ",string(other_mode)))
        box(plot_axes,'on')

        joint_max = max(joint_pdf,[],'all');
        if joint_max <= 0 || ~isfinite(joint_max)
            joint_max = 1;
        end
        clim(plot_axes,[0 joint_max])

        joint_colorbar = colorbar(plot_axes);
        joint_colorbar.Label.String = 'Probability density';
    end

    title(plot_layout,append("Joint distributions with m = ",string(mode_index)))
    xlabel(plot_layout,append("J_{",string(mode_index),"}/P"))
    ylabel(plot_layout,'J_n/P')
end

function inverse = reciprocal_occupation(occupation)
    inverse = nan(size(occupation));
    positive = isfinite(occupation) & occupation > 0;
    inverse(positive) = 1./occupation(positive);
end

function limits = tight_plot_limits(values)
    values = values(isfinite(values));
    minimum = min(values);
    maximum = max(values);

    limits = [minimum - 0.05*abs(minimum), ...
        maximum + 0.05*abs(maximum)];
    if limits(1) == limits(2)
        padding = max(0.05*abs(minimum),0.05);
        limits = limits + [-padding padding];
    end
end

function [energy_plot,occupation_plot] = ...
        connect_energy_bands(energy,occupation)
    [energy,I_order] = sort(energy(:));
    occupation = occupation(I_order);
    N_modes = numel(energy);
    break_after = false(max(N_modes - 1,0),1);

    if N_modes >= 3
        energy_span = energy(end) - energy(1);
        energy_scale = max(1,max(abs(energy)));
        tolerance = max(100*eps(energy_scale),1e-10*energy_span);
        gaps = diff(energy);
        sorted_gaps = sort(max(gaps,tolerance),'descend');
        gap_ratios = sorted_gaps(1:end-1)./sorted_gaps(2:end);
        I_cut = find(gap_ratios >= 4,1,'first');

        if ~isempty(I_cut)
            gap_threshold = sqrt(sorted_gaps(I_cut)*sorted_gaps(I_cut + 1));
            break_after = gaps > gap_threshold;
        end
    end

    N_breaks = sum(break_after);
    energy_plot = nan(N_modes + N_breaks,1);
    occupation_plot = nan(N_modes + N_breaks,1);
    jj = 1;
    for mm = 1:N_modes
        energy_plot(jj) = energy(mm);
        occupation_plot(jj) = occupation(mm);
        jj = jj + 1;
        if mm < N_modes && break_after(mm)
            jj = jj + 1;
        end
    end
end

function bar_plot = gradient_bar(ax,x,height,width,colors)
    N_colors = size(colors,1);
    y = linspace(0,height,N_colors).';
    vertices = zeros(2*N_colors,2);
    vertices(1:2:end,:) = [x - width/2 + zeros(N_colors,1),y];
    vertices(2:2:end,:) = [x + width/2 + zeros(N_colors,1),y];
    faces = [(1:2:2*N_colors - 3).',(2:2:2*N_colors - 2).', ...
        (4:2:2*N_colors).',(3:2:2*N_colors - 1).'];

    bar_plot = patch(ax,'Faces',faces,'Vertices',vertices, ...
        'FaceVertexCData',repelem(colors,2,1), ...
        'FaceColor','interp','FaceAlpha',0.92,'EdgeColor','none');
    top_color = colors(end,:) + 0.35*(1 - colors(end,:));
    plot(ax,[x - 0.42*width x + 0.42*width],[height height], ...
        'Color',top_color,'LineWidth',0.8,'HandleVisibility','off')
end

function histogram_plot = gradient_histogram( ...
        ax,edges,heights,bottom_color,top_color)
    histogram_plot = gobjects(1);
    for bb = 1:numel(heights)
        x = [edges(bb) edges(bb + 1) edges(bb + 1) edges(bb)];
        y = [0 0 heights(bb) heights(bb)];
        vertex_colors = [bottom_color;bottom_color;top_color;top_color];
        bar_plot = patch(ax,'XData',x,'YData',y, ...
            'FaceVertexCData',vertex_colors, ...
            'FaceColor','interp','EdgeColor','none');
        if bb == 1
            histogram_plot = bar_plot;
        else
            bar_plot.HandleVisibility = 'off';
        end
    end
end

function triangular_prism(ax,mode_limits,strength_limits)
    x_min = mode_limits(1);
    x_max = mode_limits(2);
    y_min = mode_limits(1);
    y_max = mode_limits(2);
    z_min = strength_limits(1);
    z_max = strength_limits(2);
    base = [
        x_min y_min z_min;
        x_max y_min z_min;
        x_max y_max z_min;
        x_min y_min z_min
    ];
    top = base;
    top(:,3) = z_max;

    hold(ax,'on')
    plot3(ax,base(:,1),base(:,2),base(:,3),'k-', ...
        'LineWidth',0.85,'HandleVisibility','off')
    plot3(ax,top(:,1),top(:,2),top(:,3),'k-', ...
        'LineWidth',0.85,'HandleVisibility','off')
    for kk = 1:3
        plot3(ax,[base(kk,1) top(kk,1)], ...
            [base(kk,2) top(kk,2)],[z_min z_max],'k-', ...
            'LineWidth',0.85,'HandleVisibility','off')
    end
    hold(ax,'off')
end
