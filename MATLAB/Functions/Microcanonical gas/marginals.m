function [x,rho] = marginals(k,x,U,P,epsilon)
    epsilon = epsilon(:).';
    x = x(x >= 0 & x < P);

    %% Isolate level of interest
    epsilon_k = epsilon(k);
    epsilon(k) = [];

    %% Slice volume
    if isscalar(U)
        U_s = U - epsilon_k*x;
        P_s = P - x;

        %% Physical slice support
        support = physical_support(U_s,P_s,epsilon);
        A = zeros(size(U_s));
        A(support) = volume( ...
            U_s(support),P_s(support),epsilon,[0 0]);

        %% Normalize
        Z = trapz(x,A);
        if Z > 0
            rho = A/Z;
        else
            rho = A;
        end
    else
        if ~isvector(U) || ~isscalar(P)
            error('marginals:BatchSize',...
                'Batched U must be a vector and P must be scalar.')
        end

        x_row = x(:).';
        U_s = U(:) - epsilon_k*x_row;
        P_s = repmat(P - x_row,numel(U),1);

        %% Physical slice support
        support = physical_support(U_s,P_s,epsilon);
        A = zeros(size(U_s));
        A(support) = volume( ...
            U_s(support),P_s(support),epsilon,[0 0]);

        %% Normalize each energy independently
        rho = A;
        for uu = 1:numel(U)
            Z = trapz(x_row,A(uu,:));
            if Z > 0
                rho(uu,:) = A(uu,:)/Z;
            end
        end
    end
end

function support = physical_support(U,P,epsilon)
    lower = min(epsilon).*P;
    upper = max(epsilon).*P;

    scale = max(abs(U),max(abs(lower),abs(upper)));
    tolerance = 32*eps(scale);

    support = isfinite(U) & isfinite(P) & P > 0 ...
        & U >= lower - tolerance ...
        & U <= upper + tolerance;
end
