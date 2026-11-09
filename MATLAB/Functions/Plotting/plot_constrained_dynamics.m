function fig = plot_constrained_dynamics(z,c,epsilon,options)
%PLOT_CONSTRAINED_DYNAMICS Animate projected occupations and relative phases.
%   FIG = PLOT_CONSTRAINED_DYNAMICS(Z,C,EPSILON) accepts the saved modal
%   amplitudes C(mode,step,realization), with at least four modes, and one
%   eigenvalue per mode. Both views share
%   propagation distance and fading tails. Occupations are red and phases
%   are blue. The displayed coordinates remain J1,J2,J3 and the phases of
%   modes 1:3 relative to the last mode M. Phases wrap
%   into [-pi,pi]; tails break at wrapping boundaries.
%
%   The occupation reference is the nonnegative polytope sum(J)=P and
%   EPSILON*J=U, projected into (J1,J2,J3), using initial ensemble-mean
%   P and linear U. All modes contribute to both constraints. Omitted
%   occupations are allowed to vary, so larger lattices generally project
%   to a 3-D volume, whose boundary is drawn translucently. Noisy members
%   need not share these values. In Kerr
%   dynamics linear U is not the conserved total energy, so trajectories
%   can leave this reference region and are never forced onto it.
%
%   Name-value options: TailLength (z units, default 50), PlaybackDuration
%   (seconds, default 30), FrameRate (default 30), AutoPlay (default true),
%   RotateOccupations (default true), and Visible ('on' or 'off'). The
%   occupation view turns once over the full z interval; the phase view
%   stays fixed. The Rotate J checkbox toggles rotation. Pause to inspect
%   or scrub; manual viewing angles are preserved on resume. Play at
%   the end replays. Closing/deleting FIG also disposes of its playback timer.
%   The Ensemble average button switches both views to one larger mean point
%   and its fading tail. Occupations use arithmetic means; phases use equal-
%   weight circular means of valid member points. If any mean phase direction
%   is undefined (zero resultant or no valid members), that point is hidden.
%   FIG.UserData exposes the reference vertices and graphics for inspection.

    arguments
        z (:,1) double {mustBeReal,mustBeFinite}
        c (:,:,:) double {mustBeFinite}
        epsilon (1,:) double {mustBeReal,mustBeFinite}
        options.TailLength (1,1) double {mustBePositive,mustBeFinite} = 50
        options.PlaybackDuration (1,1) double {mustBePositive,mustBeFinite} = 30
        options.FrameRate (1,1) double {mustBePositive,mustBeFinite} = 30
        options.AutoPlay (1,1) logical = true
        options.RotateOccupations (1,1) logical = true
        options.Visible (1,1) string {mustBeMember(options.Visible,["on","off"])} = "on"
    end
    if numel(z) < 2 || any(diff(z) <= 0)
        error('plot_constrained_dynamics:InvalidDistances', ...
            'z must contain at least two strictly increasing distances.')
    end
    M = size(c,1);
    if M < 4 || size(c,2) ~= numel(z) || isempty(c)
        error('plot_constrained_dynamics:InvalidAmplitudes', ...
            'c must have size (M,numel(z),N_init), with M >= 4.')
    end
    if numel(epsilon) ~= M
        error('plot_constrained_dynamics:InvalidSpectrum', ...
            'epsilon must contain one eigenvalue for each of the %d modes.',M)
    end

    %% Coordinates and initial reference manifold
    N_steps = numel(z);
    N_init = size(c,3);
    J = abs(c).^2;
    initialJ = mean(J(:,1,:),3);
    P = sum(initialJ);
    U = epsilon*initialJ;
    vertices = occupation_vertices(epsilon,P,U);
    occupations = permute(J(1:3,:,:),[2 1 3]); % step, coordinate, member
    phases = permute(angle(c(1:3,:,:).*conj(c(M,:,:))),[2 1 3]);
    % Phase is undefined when a participating amplitude vanishes.
    phaseFloor = 1e-12*max(J,[],[1 2]);
    validPhase = permute(all(J(1:3,:,:) > phaseFloor,1) & ...
        J(M,:,:) > phaseFloor,[2 1 3]);
    phases(repmat(~validPhase,1,3,1)) = NaN;
    averageOccupations = mean(occupations,3);
    phaseResultant = mean(exp(1i*phases),3,'omitnan');
    averagePhases = angle(phaseResultant);
    undefinedAverage = any(~isfinite(phaseResultant) | abs(phaseResultant) <= 1e-12,2);
    averagePhases(undefinedAverage,:) = NaN;

    %% Figure and shared playback controls
    ink = [0.16 0.23 0.27];
    background = [1 1 1];
    occupationColor = [0.85 0.12 0.12];
    phaseColor = [0.10 0.35 0.85];
    frame = 1;
    tailLength = options.TailLength;
    playbackTimer = [];
    clockStart = [];
    playbackStartZ = z(1);
    rotationEnabled = options.RotateOccupations;
    rotationZ = z(1);
    showAverage = false;

    fig = uifigure('Name',sprintf('Constrained %d-mode dynamics',M), ...
        'Tag','ConstrainedDynamicsFigure','Position',[80 100 1260 690], ...
        'Color',background,'Visible','off');
    fig.DeleteFcn = @dispose_timer;
    layout = uigridlayout(fig,[3 1]);
    layout.RowHeight = {'1x',42,26};
    layout.Padding = [14 8 14 8];
    layout.RowSpacing = 6;
    layout.BackgroundColor = background;
    plots = uigridlayout(layout,[1 2]);
    plots.Layout.Row = 1;
    plots.Padding = [0 0 0 0];
    plots.ColumnSpacing = 16;
    plots.BackgroundColor = background;
    axesJ = uiaxes(plots);
    axesPhase = uiaxes(plots);
    axesJ.Layout.Column = 1;
    axesPhase.Layout.Column = 2;
    for ax = [axesJ axesPhase]
        hold(ax,'on')
        grid(ax,'on')
        box(ax,'on')
        view(ax,135,24)
        pbaspect(ax,[1 1 1])
        ax.FontSize = 16;
        ax.Color = background;
        ax.XColor = ink;
        ax.YColor = ink;
        ax.ZColor = ink;
        ax.GridAlpha = 0.12;
        ax.Projection = 'perspective';
        ax.SortMethod = 'depth';
    end
    title(axesJ,'Modal occupations','FontSize',22,'Color',ink)
    subtitle(axesJ,sprintf('Projected initial constraints: P = %.4g, U = %.4g',P,U),'FontSize',15)
    xlabel(axesJ,'J_1')
    ylabel(axesJ,'J_2')
    zlabel(axesJ,'J_3')
    draw_manifold(axesJ,vertices(:,1:3));
    upper = max([reshape(permute(occupations,[1 3 2]),[],3); vertices(:,1:3)],[],1);
    upper = max(upper,max(max(upper),1)*0.05)*1.06;
    xlim(axesJ,[0 upper(1)])
    ylim(axesJ,[0 upper(2)])
    zlim(axesJ,[0 upper(3)])

    title(axesPhase,'Relative phases','FontSize',22,'Color',ink)
    subtitle(axesPhase,sprintf('Mode %d reference; angles wrap at \\pm\\pi',M),'FontSize',15)
    xlabel(axesPhase,sprintf('\\phi_1 - \\phi_{%d}',M))
    ylabel(axesPhase,sprintf('\\phi_2 - \\phi_{%d}',M))
    zlabel(axesPhase,sprintf('\\phi_3 - \\phi_{%d}',M))
    set(axesPhase,'XLim',[-pi pi],'YLim',[-pi pi],'ZLim',[-pi pi], ...
        'XTick',[-pi 0 pi],'YTick',[-pi 0 pi],'ZTick',[-pi 0 pi], ...
        'XTickLabel',{'-\pi','0','\pi'}, ...
        'YTickLabel',{'-\pi','0','\pi'}, ...
        'ZTickLabel',{'-\pi','0','\pi'});

    tailsJ = gobjects(N_init,1);
    tailsPhase = gobjects(N_init,1);
    for member = 1:N_init
        tailsJ(member) = make_tail(axesJ,occupationColor);
        tailsPhase(member) = make_tail(axesPhase,phaseColor);
    end
    headsJ = scatter3(axesJ,nan(1,N_init),nan(1,N_init),nan(1,N_init), ...
        48,occupationColor,'filled','MarkerEdgeColor',[1 1 1],'LineWidth',0.65);
    headsPhase = scatter3(axesPhase,nan(1,N_init),nan(1,N_init),nan(1,N_init), ...
        48,phaseColor,'filled','MarkerEdgeColor',[1 1 1],'LineWidth',0.65);

    controls = uigridlayout(layout,[1 9]);
    controls.Layout.Row = 2;
    controls.ColumnWidth = {78,78,20,'1x',145,44,76,100,160};
    controls.Padding = [0 3 0 3];
    controls.ColumnSpacing = 12;
    controls.BackgroundColor = background;
    playButton = uibutton(controls,'Text','Play','FontSize',14,'ButtonPushedFcn',@toggle_play);
    restartButton = uibutton(controls,'Text','Restart','FontSize',14,'ButtonPushedFcn',@restart);
    uilabel(controls,'Text','z','FontSize',16);
    distanceSlider = uislider(controls,'Limits',[z(1) z(end)],'Value',z(1), ...
        'MajorTicks',[],'MinorTicks',[], ...
        'ValueChangingFcn',@scrub,'ValueChangedFcn',@scrub);
    distanceText = uilabel(controls,'FontSize',16,'HorizontalAlignment','center');
    uilabel(controls,'Text','Tail z','FontSize',16);
    tailField = uieditfield(controls,'numeric','Value',tailLength,'FontSize',16, ...
        'Limits',[0 Inf],'LowerLimitInclusive','off','ValueChangedFcn',@change_tail);
    rotateCheckbox = uicheckbox(controls,'Text','Rotate J','FontSize',14, ...
        'Value',rotationEnabled,'ValueChangedFcn',@toggle_rotation, ...
        'Tooltip','Rotate the occupation view once over the full propagation interval.');
    averageButton = uibutton(controls,'state','Text','Ensemble average', ...
        'FontSize',14,'Value',showAverage,'ValueChangedFcn',@toggle_average, ...
        'Tooltip','Show one mean point and tail in each view. Phases use circular averages.');
    footer = uilabel(layout,'FontColor',[0.40 0.44 0.46], ...
        'Text',sprintf('%d trajectories',N_init),'FontSize',16, ...
        'HorizontalAlignment','center');
    footer.Layout.Row = 3;

    playbackTimer = timer('ExecutionMode','fixedSpacing','BusyMode','drop', ...
        'Period',max(0.001,round(1000/options.FrameRate)/1000), ...
        'Name','ConstrainedDynamicsPlayback','TimerFcn',@advance, ...
        'ErrorFcn',@timer_failed);
    fig.UserData = struct('P',P,'U',U,'ReferenceVertices',vertices,'PhaseReference',M, ...
        'InitialP',reshape(sum(J(:,1,:),1),1,[]), ...
        'InitialU',epsilon*reshape(J(:,1,:),M,[]), ...
        'Axes',[axesJ axesPhase],'Heads',[headsJ headsPhase], ...
        'OccupationTails',tailsJ,'PhaseTails',tailsPhase, ...
        'Slider',distanceSlider,'PlayButton',playButton, ...
        'RestartButton',restartButton,'TailField',tailField, ...
        'RotateCheckbox',rotateCheckbox,'AverageButton',averageButton, ...
        'Timer',playbackTimer);
    render_frame()
    fig.Visible = options.Visible;
    if options.AutoPlay
        start_playback()
    end

    function render_frame()
        tailStart = find(z >= z(frame)-tailLength,1);
        indices = tailStart:frame;
        % Fade by propagation age, independent of the sampling interval.
        alpha = 0.88*max(0,1-(z(frame)-z(indices))/tailLength).^1.7;
        if showAverage
            update_tail(tailsJ(1),averageOccupations(indices,:),alpha,false)
            update_tail(tailsPhase(1),averagePhases(indices,:),alpha,true)
            pointJ = averageOccupations(frame,:).';
            pointPhase = averagePhases(frame,:).';
        else
            for memberIndex = 1:N_init
                update_tail(tailsJ(memberIndex),occupations(indices,:,memberIndex),alpha,false)
                update_tail(tailsPhase(memberIndex),phases(indices,:,memberIndex),alpha,true)
            end
            pointJ = reshape(occupations(frame,:,:),3,N_init);
            pointPhase = reshape(phases(frame,:,:),3,N_init);
        end
        set(headsJ,'XData',pointJ(1,:),'YData',pointJ(2,:),'ZData',pointJ(3,:))
        set(headsPhase,'XData',pointPhase(1,:),'YData',pointPhase(2,:),'ZData',pointPhase(3,:))
        distanceSlider.Value = z(frame);
        distanceText.Text = sprintf('%.3g / %.3g',z(frame),z(end));
        % Advance with propagation distance, not wall time or frame count.
        % Starting from the current view preserves manual camera adjustments.
        if rotationEnabled && z(frame) ~= rotationZ
            [azimuth,elevation] = view(axesJ);
            deltaAzimuth = 360*(z(frame)-rotationZ)/(z(end)-z(1));
            view(axesJ,azimuth+deltaAzimuth,elevation)
        end
        rotationZ = z(frame);
        setappdata(fig,'FrameIndex',frame)
    end

    function start_playback()
        playbackStartZ = z(frame);
        clockStart = tic;
        playButton.Text = 'Pause';
        start(playbackTimer)
    end

    function pause_playback()
        stop(playbackTimer)
        playButton.Text = 'Play';
    end

    function toggle_play(~,~)
        if strcmp(playbackTimer.Running,'on')
            pause_playback()
        else
            if frame == N_steps
                frame = 1;
                render_frame()
            end
            start_playback()
        end
    end

    function advance(~,~)
        if ~isgraphics(fig)
            return
        end
        currentZ = playbackStartZ + toc(clockStart)* ...
            (z(end)-z(1))/options.PlaybackDuration;
        nextFrame = find(z <= min(currentZ,z(end)),1,'last');
        if nextFrame ~= frame
            frame = nextFrame;
            render_frame()
        end
        if frame == N_steps
            pause_playback()
        end
    end

    function scrub(~,event)
        pause_playback()
        [~,frame] = min(abs(z-event.Value));
        render_frame()
    end

    function restart(~,~)
        pause_playback()
        frame = 1;
        render_frame()
        start_playback()
    end

    function change_tail(~,~)
        tailLength = tailField.Value;
        render_frame()
    end

    function toggle_rotation(~,~)
        rotationEnabled = rotateCheckbox.Value;
    end

    function toggle_average(~,~)
        showAverage = averageButton.Value;
        if showAverage
            set([tailsJ(2:end);tailsPhase(2:end)],'Visible','off')
            set([tailsJ(1);tailsPhase(1)],'LineWidth',2.4)
            set([headsJ headsPhase],'SizeData',180)
            footer.Text = sprintf('Ensemble average of %d trajectories',N_init);
        else
            set([tailsJ;tailsPhase],'Visible','on','LineWidth',1.7)
            set([headsJ headsPhase],'SizeData',48)
            footer.Text = sprintf('%d trajectories',N_init);
        end
        render_frame()
    end

    function timer_failed(~,~)
        if isgraphics(fig)
            pause_playback()
        end
    end

    function dispose_timer(~,~)
        if ~isempty(playbackTimer) && isvalid(playbackTimer)
            stop(playbackTimer)
            delete(playbackTimer)
        end
    end
end

function handle = make_tail(ax,color)
    handle = patch(ax,'Vertices',nan(3,3),'Faces',[1 2 3], ...
        'FaceColor','none','EdgeColor',color,'EdgeAlpha','interp', ...
        'FaceVertexAlphaData',zeros(3,1),'AlphaDataMapping','none', ...
        'LineWidth',1.7,'HitTest','off','PickableParts','none');
end

function update_tail(handle,points,alpha,wrapped)
    % Each open edge ends at a NaN vertex so patch never closes the tail.
    N = size(points,1);
    if N < 2
        set(handle,'Faces',[1 2 3],'Vertices',nan(3,3), ...
            'FaceVertexAlphaData',zeros(3,1))
        return
    end
    edges = [(1:N-1).' (2:N).'];
    keep = all(isfinite(points(1:end-1,:)),2) & all(isfinite(points(2:end,:)),2);
    if wrapped
        keep = keep & all(abs(diff(points,1,1)) <= pi,2);
    end
    edges = edges(keep,:);
    if isempty(edges)
        faces = [N+1 N+1 N+1];
    else
        faces = [edges repmat(N+1,size(edges,1),1)];
    end
    set(handle,'Faces',faces,'Vertices',[points;nan(1,3)], ...
        'FaceVertexAlphaData',[alpha(:);0])
end

function vertices = occupation_vertices(epsilon,P,U)
    % An extreme point under two linear equalities has <=2 occupied modes.
    M = numel(epsilon);
    if P == 0
        vertices = zeros(1,M);
        return
    end
    span = max(epsilon)-min(epsilon);
    if span == 0
        vertices = P*eye(M); % energy constraint is redundant
        return
    end
    e = (epsilon-min(epsilon))/span;
    target = (U/P-min(epsilon))/span;
    tolerance = 64*eps;
    onTarget = find(abs(e-target) <= tolerance);
    below = find(e < target-tolerance);
    above = find(e > target+tolerance);
    vertices = zeros(numel(onTarget)+numel(below)*numel(above),M);
    row = 0;
    for mode = onTarget
        row = row+1;
        vertices(row,mode) = P;
    end
    for low = below
        for high = above
            row = row+1;
            fraction = (e(high)-target)/(e(high)-e(low));
            vertices(row,low) = P*fraction;
            vertices(row,high) = P-vertices(row,low);
        end
    end
end

function draw_manifold(ax,xyz)
    % Projection commutes with convex hull. Distinct full-space vertices
    % can coincide when their occupied modes are all outside the plot.
    xyz = unique(xyz,'rows','stable');
    color = [0.24 0.55 0.59];
    centered = xyz-mean(xyz,1);
    [~,S,basis] = svd(centered,0);
    dimension = sum(diag(S) > max(size(centered))*eps(max(norm(centered),1)));
    if dimension == 0
        plot3(ax,xyz(1,1),xyz(1,2),xyz(1,3),'o','Color',color)
    elseif dimension == 1
        [~,order] = sort(centered*basis(:,1));
        plot3(ax,xyz(order,1),xyz(order,2),xyz(order,3),'-','Color',color,'LineWidth',1.4)
    else
        if dimension == 2
            plane = centered*basis(:,1:2);
            order = convhull(plane(:,1),plane(:,2));
            faces = order(1:end-1).';
        else
            faces = convhull(xyz); % boundary of the projected 3-D region
        end
        patch(ax,'Vertices',xyz,'Faces',faces,'FaceColor',color, ...
            'FaceAlpha',0.13,'EdgeColor',color,'EdgeAlpha',0.5,'LineWidth',1.3, ...
            'HitTest','off','PickableParts','none');
    end
end
