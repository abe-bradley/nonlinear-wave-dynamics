function rho = pearson_correlations(U,P,epsilon)
    epsilon = epsilon(:).';

    M = numel(epsilon);
    [~,f,d] = degeneracy(epsilon);
    N_levels = numel(f);

    %% One representative per degenerate level
    mu = zeros(1,N_levels);
    second_moment = zeros(1,N_levels);

    for aa = 1:N_levels
        alpha = zeros(1,M);
        ii = f(aa);

        alpha(ii) = 1;
        mu(aa) = mixed_moments(U,P,epsilon,alpha);

        alpha(ii) = 2;
        second_moment(aa) = mixed_moments(U,P,epsilon,alpha);
    end

    variance = second_moment - mu.^2;
    sigma = sqrt(variance);

    %% Equal-energy blocks
    rho = zeros(M);
    for aa = 1:N_levels
        indices = d{f(aa)};

        % For distinct modes in the same level, the augmented pole orders
        % are identical to the self moment, but without its factor of 2!.
        rho_same = (second_moment(aa)/2 - mu(aa)^2)/variance(aa);
        block = rho_same*ones(numel(indices),numel(indices));
        block(1:numel(indices)+1:end) = 1;
        rho(indices,indices) = block;
    end

    %% Cross-level blocks
    for aa = 1:N_levels
        for bb = aa+1:N_levels
            alpha = zeros(1,M);
            alpha(f(aa)) = 1;
            alpha(f(bb)) = 1;

            cross_moment = mixed_moments(U,P,epsilon,alpha);
            covariance = cross_moment - mu(aa)*mu(bb);
            rho_ab = covariance/(sigma(aa)*sigma(bb));

            indices_a = d{f(aa)};
            indices_b = d{f(bb)};
            rho(indices_a,indices_b) = rho_ab;
            rho(indices_b,indices_a) = rho_ab;
        end
    end
end
