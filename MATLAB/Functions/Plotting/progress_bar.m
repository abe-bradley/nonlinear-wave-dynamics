function fig = progress_bar(z,zf,fig)
    %PROGRESS_BAR Floating progress bar with an internal completion timer.
    %
    % fig = progress_bar(z,zf,fig) initializes or updates a progress bar.
    % The first call should pass [] for fig. Later calls should pass the
    % returned figure handle back into the function.
    %
    % Inputs:
    % z   - Current progress coordinate. It can be any monotonically
    %       increasing distance-like quantity.
    % zf  - Final progress coordinate, in the same units as z.
    % fig - Existing progress-bar figure handle. Pass [] to initialize.
    %
    % Output:
    % fig - Progress-bar figure handle. Reuse this handle on later calls.
    %
    % When z >= zf, the bar fills, the total elapsed time is printed, and
    % the figure closes automatically.

    %EXAMPLE USAGE FOR CONTINUOUS ITERATIONS
    % progress_fig = progress_bar(0,Zf,[]);
    % for nn = 1:length(Z)
    %     [code block]
    %     progress_fig = progress_bar(Z(nn-1),Zf,progress_fig);
    % end
    % progress_fig = progress_bar(Zf,Zf,progress_fig);

    %EXAMPLE USAGE FOR DISCRETE ITERATIONS
    % progress_fig = progress_bar(0,N,[]);
    % for nn = 1:N
    %     [code block]
    %     progress_fig = progress_bar(nn,N,progress_fig);
    % end


    narginchk(3,3)

    new_fig = isempty(fig) || ~isvalid(fig) || ~isstruct(fig.UserData) || ~isfield(fig.UserData,'progress_bar');
    if new_fig
        fig = initialize_progress_bar(z);
    else
        % Record the time and distance covered since the previous update.
        handles = fig.UserData.progress_bar;
        dt = toc(handles.tic);
        dz = z - handles.z;

        if dz > 0 && dt > 0
            handles.N_updates = handles.N_updates + 1;
            handles.time = [handles.time dt];
            handles.distance = [handles.distance dz];
        end

        handles.tic = tic;
        handles.z = z;
        fig.UserData.progress_bar = handles;
    end

    handles = fig.UserData.progress_bar;
    prog = min(1,max(0,z/zf));
    done = z >= zf;

    if done
        time_string = format_duration(0);
    elseif isempty(handles.distance)
        time_string = '--h --m --s';
    else
        % Recent updates receive more weight in the remaining-time estimate.
        weights = 1:numel(handles.time);
        time_block = sum(weights.*handles.time)/sum(weights);
        distance_block = sum(weights.*handles.distance)/sum(weights);
        remaining_time = time_block*max(0,zf - z)/distance_block;
        time_string = format_duration(remaining_time);
    end

    handles.title.String = ['Estimated time remaining: ' time_string];
    if prog == 0
        handles.fill.Visible = 'off';
    else
        handles.fill.Visible = 'on';
        handles.fill.Position = [handles.bar(1) handles.bar(2) handles.bar(3)*prog handles.bar(4)];
    end

    drawnow limitrate

    if new_fig
        % Start timing after the initial graphics draw.
        handles = fig.UserData.progress_bar;
        handles.tic = tic;
        fig.UserData.progress_bar = handles;
    end

    handles = fig.UserData.progress_bar;
    if done && ~handles.done
        fprintf('Computation completed in %s.\n',format_duration(sum(handles.time)))
        handles.done = true;
        fig.UserData.progress_bar = handles;
        close(fig)
        fig = [];
    end
end

function fig = initialize_progress_bar(z)
    fig = figure('Name','Propagation progress','NumberTitle','off', ...
        'MenuBar','none','ToolBar','none','Color','w','Position',[300 300 760 150]);

    handles.bar = [0.08 0.50 0.84 0.18];
    handles.fill = annotation(fig,'rectangle',[handles.bar(1) handles.bar(2) eps handles.bar(4)], ...
        'FaceColor',[0.28 0.68 0.18],'EdgeColor','none','Visible','off');
    handles.frame = annotation(fig,'rectangle',handles.bar, ...
        'FaceColor','none','Color','k','LineWidth',1.5);

    handles.title = uicontrol(fig,'Style','text','Units','normalized', ...
        'Position',[0.08 0.74 0.84 0.18],'BackgroundColor','w', ...
        'FontSize',13,'String','Estimated time remaining: --h --m --s');

    ticks = [0 0.25 0.5 0.75 1];
    labels = {'0%','25%','50%','75%','100%'};
    text_width = 0.14;
    for kk = 1:numel(ticks)
        x_tick = handles.bar(1) + handles.bar(3)*ticks(kk);
        annotation(fig,'line',[x_tick x_tick],[0.34 0.43],'Color','k','LineWidth',1);
        uicontrol(fig,'Style','text','Units','normalized', ...
            'Position',[x_tick - text_width/2 0.12 text_width 0.18], ...
            'BackgroundColor','w','FontSize',13, ...
            'HorizontalAlignment','center','String',labels{kk});
    end

    handles.tic = tic;
    handles.z = z;
    handles.N_updates = 0;
    handles.time = [];
    handles.distance = [];
    handles.done = false;

    fig.UserData.progress_bar = handles;
end

function str = format_duration(t)
    t = max(0,round(t));
    hours = floor(t/3600);
    minutes = floor(mod(t,3600)/60);
    seconds = mod(t,60);

    str = sprintf('%02dh %02dm %02ds',hours,minutes,seconds);
end
