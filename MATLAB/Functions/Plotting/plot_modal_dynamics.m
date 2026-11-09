function plot_modal_dynamics(N_end,z,c,sysk)
    %% Inputs
    narginchk(4,4)

    [c,c_hierarchical] = system_hierarchy(c);
    [sysk,sysk_hierarchical] = system_hierarchy(sysk);
    if c_hierarchical ~= sysk_hierarchical
        error('c and sysk must use the same system hierarchy.')
    end
    if ~isscalar(c) || ~isscalar(sysk) || ~isstruct(sysk.full) ...
            || ~isscalar(sysk.full)
        error('c and sysk must be scalar system containers.')
    end

    z = z(:);
    N_steps = numel(z);
    c.full = modal_dynamics_array(c.full,N_steps,'c.full');

    M = size(c.full,1);
    N_init = size(c.full,3);
    if ~isscalar(N_end) || N_end < 1 || N_end > N_steps || N_end ~= round(N_end)
        error('N_end must be an integer between 1 and the number of propagation steps.')
    end

    %% Subsystems
    B = cell(0,1);
    if sysk_hierarchical
        if ~isfield(sysk,'B') || ~iscell(sysk.B) || isempty(sysk.B) ...
                || ~isstruct(sysk.sub) || ~isscalar(sysk.sub)
            error('A hierarchical sysk must contain packed sub and nonempty B.')
        end
        B = sysk.B;
    end

    N_sub = numel(B);
    if N_sub > 0
        c.sub = modal_dynamics_array(c.sub,N_steps,'c.sub');
        if size(c.sub,1) ~= M || size(c.sub,3) ~= N_init
            error('c.sub must have the same dimensions as c.full.')
        end

        all_indices = zeros(1,0);
        for ss = 1:N_sub
            b = B{ss};
            if isempty(b) || ~isvector(b) || any(b ~= round(b)) ...
                    || any(b < 1) || any(b > M)
                error('B{%d} contains invalid modal indices.',ss)
            end
            B{ss} = b(:).';
            all_indices = [all_indices B{ss}]; %#ok<AGROW>
        end
        if ~isequal(all_indices,1:M)
            error('B must partition the packed subsystem modes in order.')
        end
    end

    if N_sub <= 1
        N_sub_plot = 0;
    else
        N_sub_plot = N_sub;
    end

    %% Figure
    viewer = uifigure( ...
        'Name','Modal dynamics', ...
        'Tag','ModalDynamicsFigure', ...
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

    toolbar_grid = uigridlayout(toolbar,[1 5]);
    toolbar_grid.ColumnWidth = {70,220,55,260,'1x'};
    toolbar_grid.RowHeight = {'1x'};
    toolbar_grid.Padding = [12 7 12 7];
    toolbar_grid.ColumnSpacing = 8;

    system_label = uilabel(toolbar_grid, ...
        'Text','System', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center');
    system_label.Layout.Column = 1;

    system_names = [{'Full system'}, ...
        arrayfun(@(ss) sprintf('Subsystem %d',ss),1:N_sub_plot,'UniformOutput',false)];
    system_dropdown = uidropdown(toolbar_grid, ...
        'Items',system_names, ...
        'ItemsData',1:numel(system_names), ...
        'Value',1, ...
        'FontSize',16, ...
        'Tag','ModalDynamicsSystemDropdown');
    system_dropdown.Layout.Column = 2;

    plot_label = uilabel(toolbar_grid, ...
        'Text','Plot', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','center');
    plot_label.Layout.Column = 3;

    plot_dropdown = uidropdown(toolbar_grid, ...
        'Items',{'Modal evolution','Entropy evolution', ...
        'Three-mode trajectories','Complex planes','Mutual coherence'}, ...
        'ItemsData',{'evolution','entropy','trajectories','complex','coherence'}, ...
        'Value','evolution', ...
        'FontSize',16, ...
        'Tag','ModalDynamicsPlotDropdown');
    plot_dropdown.Layout.Column = 4;

    %% Plot panel
    plot_panel = uipanel(viewer_grid, ...
        'BackgroundColor',[1 1 1], ...
        'BorderType','none', ...
        'Tag','ModalDynamicsPlotPanel');
    plot_panel.Layout.Row = 2;
    plot_panel.Layout.Column = 1;

    %% Callbacks
    system_dropdown.ValueChangedFcn = @selection_changed;
    plot_dropdown.ValueChangedFcn = @selection_changed;

    viewer.Visible = 'on';
    drawnow
    redraw_viewer()

    function selection_changed(~,~)
        redraw_viewer()
    end

    function redraw_viewer()
        if ~isvalid(viewer) || ~isvalid(plot_panel)
            return
        end

        system_dropdown.Enable = 'off';
        plot_dropdown.Enable = 'off';
        viewer.Pointer = 'watch';
        drawnow
        cleanup = onCleanup(@() modal_dynamics_restore_controls( ...
            viewer,system_dropdown,plot_dropdown));

        delete(plot_panel.Children)

        system_index = system_dropdown.Value;
        if system_index == 1
            c_plot = c.full;
        else
            c_plot = c.sub(B{system_index - 1},:,:);
        end

        switch string(plot_dropdown.Value)
            case "evolution"
                plot_modal_evolution(plot_panel,z,c_plot)
            case "entropy"
                plot_entropy_evolution(plot_panel,z,c_plot)
            case "trajectories"
                plot_three_mode_trajectories(plot_panel,z,c_plot)
            case "complex"
                plot_complex_planes(plot_panel,c_plot,N_end)
            case "coherence"
                plot_mutual_coherence(plot_panel,c_plot,N_end)
        end

        drawnow
        clear cleanup
    end
end

function c = modal_dynamics_array(c,N_steps,name)
    if ~isnumeric(c) || isempty(c) || size(c,2) ~= N_steps || ndims(c) > 3
        error('%s must have size (N_modes,N_steps,N_init).',name)
    end
    if ismatrix(c)
        c = reshape(c,size(c,1),N_steps,1);
    end
end

function modal_dynamics_restore_controls(viewer,system_dropdown,plot_dropdown)
    if isvalid(viewer)
        viewer.Pointer = 'arrow';
    end
    if isvalid(system_dropdown)
        system_dropdown.Enable = 'on';
    end
    if isvalid(plot_dropdown)
        plot_dropdown.Enable = 'on';
    end
end

function modal_dynamics_message(parent,message)
    message_grid = uigridlayout(parent,[1 1]);
    message_grid.Padding = [0 0 0 0];
    message_label = uilabel(message_grid, ...
        'Text',message, ...
        'FontSize',20, ...
        'FontColor',[0.45 0.48 0.52], ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','center');
    message_label.Layout.Row = 1;
    message_label.Layout.Column = 1;
end

function [J,j_ave,P0,mode_colors] = modal_dynamics_data(c)
    M = size(c,1);
    J = abs(c).^2;
    J_ave = mean(J,3);
    P0 = mean(sum(J(:,1,:),1),3);
    if P0 == 0
        P0 = 1;
    end
    j_ave = J_ave/P0;

    mode_colors = ot_colormap("guppy");
    I_colors = round(linspace(1,size(mode_colors,1),M));
    mode_colors = flip(mode_colors(I_colors,:),1);
end

function plot_modal_evolution(parent,z,c)
    %% Modal occupations
    M = size(c,1);
    N_steps = size(c,2);
    [~,j_ave,~,mode_colors] = modal_dynamics_data(c);

    N_z_plot = min(N_steps,1000);
    I_plot = unique(round(linspace(1,N_steps,N_z_plot)));
    z_plot = z(I_plot);
    j_plot = j_ave(:,I_plot);

    %% Gaussian ridges
    x = linspace(-0.5,0.5,100);
    width = 0.1;
    gaussian = exp(-(x/width).^2);

    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    hold(plot_axes,'on')

    for mm = 1:M
        ridge = j_plot(mm,:).'*gaussian;

        fade = reshape(gaussian.^1.3,1,[],1);
        mode_color = reshape(mode_colors(mm,:),1,1,3);
        ridge_color = 1 - fade.*(1 - mode_color);
        ridge_color = repmat(ridge_color,numel(I_plot),1,1);

        surf(plot_axes,x + mm,z_plot,ridge,ridge_color, ...
            'FaceColor','interp','EdgeColor','none')
        plot3(plot_axes,mm*ones(size(z_plot)),z_plot,j_plot(mm,:), ...
            'Color',[0.0 0.0 0.1],'LineWidth',0.5)
    end

    hold(plot_axes,'off')
    xlabel(plot_axes,'Mode')
    ylabel(plot_axes,'z')
    zlabel(plot_axes,'J_m/P_0')
    title(plot_axes,'Ensemble-averaged modal evolution')
    xlim(plot_axes,[0.5 M + 0.5])

    occupation_max = max(j_plot,[],'all');
    if ~isfinite(occupation_max) || occupation_max <= 0
        occupation_max = 1;
    end
    zlim(plot_axes,[0 max(1,1.05*occupation_max)])
    grid(plot_axes,'on')
    box(plot_axes,'on')
    pbaspect(plot_axes,[1.0 1.0 0.8])
    view(plot_axes,[-135 15])
    set(plot_axes,'XDir','reverse')
end

function plot_entropy_evolution(parent,z,c)
    %% Entropy
    [~,j_ave] = modal_dynamics_data(c);
    S = sum(log(j_ave),1).';
    S(~isfinite(S)) = nan;

    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    plot(plot_axes,z,S,'k','LineWidth',2)

    xlabel(plot_axes,'z')
    ylabel(plot_axes,'S')
    title(plot_axes,'Entropy evolution')
    grid(plot_axes,'on')
    box(plot_axes,'on')
end

function plot_three_mode_trajectories(parent,z,c)
    %% Dimensions
    M = size(c,1);
    if M < 3
        modal_dynamics_message(parent, ...
            'Three-mode trajectories require at least three modes.')
        return
    end

    N_steps = size(c,2);
    N_init = size(c,3);
    [J,~,P0] = modal_dynamics_data(c);

    J1 = reshape(J(1,:,:),N_steps,N_init)/P0;
    J2 = reshape(J(2,:,:),N_steps,N_init)/P0;
    J3 = reshape(J(3,:,:),N_steps,N_init)/P0;

    %% Trajectories
    plot_layout = tiledlayout(parent,1,1, ...
        'TileSpacing','compact','Padding','compact');
    plot_axes = nexttile(plot_layout);
    hold(plot_axes,'on')

    trajectory_colors = flip(ot_colormap("thermal-2"),1);
    colormap(plot_axes,trajectory_colors)

    for kk = 1:N_init
        patch(plot_axes,[J1(:,kk).' nan],[J2(:,kk).' nan], ...
            [J3(:,kk).' nan],[z(:).' nan], ...
            'FaceColor','none','EdgeColor','interp','LineWidth',0.75)
    end

    hold(plot_axes,'off')
    xlabel(plot_axes,'J_1/P_0')
    ylabel(plot_axes,'J_2/P_0')
    zlabel(plot_axes,'J_3/P_0')
    title(plot_axes,'Three-mode occupation-space trajectories')
    grid(plot_axes,'on')
    box(plot_axes,'on')
    view(plot_axes,[135 25])

    occupation_scale = max([J1;J2;J3],[],'all');
    if ~isfinite(occupation_scale) || occupation_scale <= 0
        occupation_scale = 1;
    end
    xlim(plot_axes,[0 max(0.05*occupation_scale,1.05*max(J1,[],'all'))])
    ylim(plot_axes,[0 max(0.05*occupation_scale,1.05*max(J2,[],'all'))])
    zlim(plot_axes,[0 max(0.05*occupation_scale,1.05*max(J3,[],'all'))])

    if z(end) > z(1)
        clim(plot_axes,[z(1) z(end)])
    end
    trajectory_colorbar = colorbar(plot_axes,'southoutside');
    trajectory_colorbar.Label.String = 'z';
end

function plot_complex_planes(parent,c,N_end)
    %% Modal occupations
    M = size(c,1);
    N_steps = size(c,2);
    N_init = size(c,3);
    I_end = N_steps - N_end + 1:N_steps;
    [~,~,P0,mode_colors] = modal_dynamics_data(c);

    c_initial = reshape(c(:,1,:),M,N_init);
    c_late = reshape(c(:,I_end,:),M,[]);

    c_initial = abs(c_initial).^2.*exp(1i*angle(c_initial))/P0;
    c_late = abs(c_late).^2.*exp(1i*angle(c_late))/P0;

    J_initial = mean(abs(c(:,1,:)).^2,3)/P0;
    J_late = mean(abs(c(:,I_end,:)).^2,[2 3])/P0;

    r_max = max(abs([c_initial(:);c_late(:)]));
    if ~isfinite(r_max) || r_max == 0
        r_max = 1;
    end
    r_limit = 1.1*max(1,r_max);

    %% Complex planes
    theta = linspace(0,2*pi,200);
    plot_grid = uigridlayout(parent,[4 1]);
    plot_grid.RowHeight = {24,'1x',24,'1x'};
    plot_grid.ColumnWidth = {'1x'};
    plot_grid.Padding = [4 4 4 4];
    plot_grid.RowSpacing = 2;

    initial_label = uilabel(plot_grid, ...
        'Text','Initial occupations', ...
        'FontSize',16, ...
        'HorizontalAlignment','center');
    initial_label.Layout.Row = 1;

    initial_grid = uigridlayout(plot_grid,[1 M]);
    initial_grid.Layout.Row = 2;
    initial_grid.ColumnWidth = repmat({'1x'},1,M);
    initial_grid.Padding = [0 0 0 0];
    initial_grid.ColumnSpacing = 2;

    late_label = uilabel(plot_grid, ...
        'Text','Late occupations', ...
        'FontSize',16, ...
        'HorizontalAlignment','center');
    late_label.Layout.Row = 3;

    late_grid = uigridlayout(plot_grid,[1 M]);
    late_grid.Layout.Row = 4;
    late_grid.ColumnWidth = repmat({'1x'},1,M);
    late_grid.Padding = [0 0 0 0];
    late_grid.ColumnSpacing = 2;

    for jj = 1:2
        for mm = 1:M
            if jj == 1
                row_grid = initial_grid;
            else
                row_grid = late_grid;
            end

            axes_panel = uipanel(row_grid, ...
                'BackgroundColor',[1 1 1], ...
                'BorderType','none');
            axes_panel.Layout.Column = mm;
            plot_axes = polaraxes(axes_panel);
            plot_axes.Units = 'normalized';
            plot_axes.Position = [0.06 0.06 0.88 0.88];
            hold(plot_axes,'on')

            if jj == 1
                c_plot = c_initial(mm,:);
                J_plot = J_initial(mm);
            else
                c_plot = c_late(mm,:);
                J_plot = J_late(mm);
            end

            polarplot(plot_axes,c_plot,'.','Color',mode_colors(mm,:), ...
                'MarkerSize',8)
            polarplot(plot_axes,theta,J_plot*ones(size(theta)), ...
                'Color',[1.0 0.0 0.4],'LineWidth',2)

            hold(plot_axes,'off')
            rlim(plot_axes,[0 r_limit])
            plot_axes.ThetaZeroLocation = 'right';
            plot_axes.ThetaDir = 'clockwise';

            if jj == 1 && mm == 1
                plot_axes.ThetaTick = 0:90:270;
                plot_axes.ThetaTickLabel = {'0','90','180',''};
                plot_axes.RAxisLocation = 135;
                plot_axes.RTick = [0 0.5 1];
                plot_axes.RTickLabel = {'0.0','0.5','1.0'};
                plot_axes.ThetaAxis.Label.String = 'arg c_m';
                plot_axes.ThetaAxis.Label.Position = [225 1.2*r_limit 0];
                plot_axes.ThetaAxis.Label.Rotation = 0;
                plot_axes.RAxis.Label.String = '|c_m|^2/P_0';
                plot_axes.RAxis.Label.Position = [135 1.2*r_limit 0];
                plot_axes.RAxis.Label.Rotation = 0;
            else
                plot_axes.RTickLabel = {};
                plot_axes.ThetaTickLabel = {};
            end

            title(plot_axes,append("m = ",string(mm)))
        end
    end
end

function plot_mutual_coherence(parent,c,N_end)
    %% Coherence matrices
    M = size(c,1);
    N_steps = size(c,2);
    N_init = size(c,3);
    I_end = N_steps - N_end + 1:N_steps;
    [~,~,P0] = modal_dynamics_data(c);

    c_initial = reshape(c(:,1,:),M,N_init);
    c_late = reshape(c(:,I_end,:),M,[]);

    Gamma_initial = (c_initial*c_initial')/(N_init*P0);
    Gamma_late = (c_late*c_late')/(size(c_late,2)*P0);
    Gamma = {Gamma_initial,Gamma_late};

    coherence_titles = {'Initial mutual coherence','Late mutual coherence'};
    coherence_max = max(abs([Gamma_initial(:);Gamma_late(:)]));
    if ~isfinite(coherence_max) || coherence_max == 0
        coherence_max = 1;
    end

    %% Initial and late coherence
    plot_layout = tiledlayout(parent,1,2, ...
        'TileSpacing','compact','Padding','compact');
    coherence_colors = sqrt(ot_colormap('Purples'));

    for jj = 1:2
        plot_axes = nexttile(plot_layout);
        coherence_bars = bar3(plot_axes,abs(Gamma{jj}));
        colormap(plot_axes,coherence_colors)

        for kk = 1:numel(coherence_bars)
            coherence_bars(kk).CData = coherence_bars(kk).ZData;
            coherence_bars(kk).FaceColor = 'interp';
        end

        xlabel(plot_axes,'n')
        ylabel(plot_axes,'m')
        zlabel(plot_axes,'|\Gamma_{mn}|/P_0')
        title(plot_axes,coherence_titles{jj})
        xlim(plot_axes,[0.5 M + 0.5])
        ylim(plot_axes,[0.5 M + 0.5])
        zlim(plot_axes,[0 1.1*coherence_max])
        clim(plot_axes,[0 coherence_max])
        grid(plot_axes,'on')
        view(plot_axes,[-35 20])
    end
end
