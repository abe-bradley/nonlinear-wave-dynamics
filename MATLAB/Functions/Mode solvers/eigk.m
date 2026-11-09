function sysk = eigk(plotting,H,B)
    % the 'k' is for Kerr!

    %% Inputs
    narginchk(2,3)

    if ~islogical(plotting) || ~isscalar(plotting)
        error('eigk:InvalidPlottingOption','Plotting must be true or false.')
    end
    if ~isnumeric(H) || isempty(H) || ~ismatrix(H) || size(H,1) ~= size(H,2)
        error('eigk:InvalidHamiltonian','H must be a nonempty square numeric matrix.')
    end

    %% Kerr properties
    if nargin < 3
        sysk = eigk_system(H);
        sysk_plot = sysk;
    else
        if ~iscell(B)
            error('eigk:InvalidSubsystems','B must be a cell array of subsystem indices.')
        end

        N = size(H,1);
        all_indices = zeros(1,0);
        for ss = 1:numel(B)
            b = B{ss};
            if ~isnumeric(b) || ~isreal(b) || isempty(b) || ~isvector(b) ...
                    || any(b ~= round(b)) ...
                    || any(b < 1) || any(b > N)
                error('eigk:InvalidSubsystems', ...
                    'B{%d} contains invalid site indices.',ss)
            end
            B{ss} = b(:).';
            all_indices = [all_indices B{ss}]; %#ok<AGROW>
        end
        if ~isequal(all_indices,1:N)
            error('eigk:InvalidSubsystems', ...
                'B must partition the packed sites in order.')
        end

        sysk.full = eigk_system(H);
        sysk.sub = eigk_subsystems(H,B);
        sysk.B = B;
        sysk_plot = sysk.full;
    end

    %% Plot Kerr interactions
    if plotting
        figure('Name','eigk K','Color','w')
        axes_handle = axes;
        imagesc(axes_handle,real(sysk_plot.K))
        title(axes_handle,'K')
        axis(axes_handle,'image')
        set(axes_handle,'YDir','normal')
        xlabel(axes_handle,'Mode index n')
        ylabel(axes_handle,'Mode index m')
        colorbar(axes_handle)
    end
end

function sysk = eigk_system(H)
    %% Eigensystem
    [psi,epsilon] = eig(H,'vector');
    [~,I_order] = sortrows([real(epsilon),imag(epsilon)],[1 2]);
    epsilon = epsilon(I_order).';
    psi = psi(:,I_order);

    %% Gram matrix
    G = psi'*psi;

    %% Kerr interactions
    N = length(epsilon);
    A = reshape(conj(psi),N,N,1).*reshape(psi,N,1,N);
    A = reshape(A,N,N^2);
    W = A.'*A;

    K = zeros(N,N);
    I_W = reshape(1:N^2,N,N);
    for mm = 1:N
        for nn = 1:N
            del = mm == nn;
            K(mm,nn) = W(I_W(mm,mm),I_W(nn,nn)) ...
                + W(I_W(mm,nn),I_W(nn,mm)) ...
                - del*W(I_W(mm,mm),I_W(mm,mm));
        end
    end
    K = K/2;

    %% Output
    sysk = struct;
    sysk.epsilon = epsilon;
    sysk.psi = psi;
    sysk.G = G;
    sysk.W = W;
    sysk.K = K;
end

function sysk = eigk_subsystems(H,B)
    %% Packed subsystem properties
    N = size(H,1);
    epsilon = zeros(1,N);
    psi = zeros(N,N);
    G = zeros(N,N);
    K = zeros(N,N);
    W = sparse(N^2,N^2);

    for ss = 1:numel(B)
        b = B{ss};
        N_sub = numel(b);
        sysk_local = eigk_system(H(b,b));

        epsilon(b) = sysk_local.epsilon;
        psi(b,b) = sysk_local.psi;
        G(b,b) = sysk_local.G;
        K(b,b) = sysk_local.K;

        S = sparse(b,1:N_sub,1,N,N_sub);
        S_pair = kron(S,S);
        W = W + S_pair*sparse(sysk_local.W)*S_pair.';
    end

    %% Output
    sysk = struct;
    sysk.epsilon = epsilon;
    sysk.psi = psi;
    sysk.G = G;
    sysk.W = W;
    sysk.K = K;
end
