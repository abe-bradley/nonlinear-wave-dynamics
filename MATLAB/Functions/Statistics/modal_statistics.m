function statistics = modal_statistics(J)
%MODAL_STATISTICS Statistics of time-dependent modal occupations.
% J has dimensions modes-by-steps-by-initial-conditions.

    %% Dimensions
    if ~isnumeric(J) || isempty(J) || ~isreal(J) ...
            || any(~isfinite(J),'all') || any(J < 0,'all') || ndims(J) > 3
        error('J must be a finite, real, nonnegative numeric array.')
    end

    [N_modes,N_steps,N_init] = size(J);
    N_membs = N_steps*N_init;
    J = reshape(double(J),N_modes,N_membs);

    statistics.N_modes = N_modes;
    statistics.N_steps = N_steps;
    statistics.N_init = N_init;
    statistics.N_membs = N_membs;
    statistics.modes = (1:N_modes).';

    %% First moments
    first_moment = mean(J,2);

    statistics.first_moment = first_moment;

    %% Second moments
    second_moment = mean(J.^2,2);
    cross_moment = (J*J.')/N_membs;
    normalized_cross_moment = cross_moment./(first_moment*first_moment.');

    statistics.second_moment = second_moment;
    statistics.cross_moment = cross_moment;
    statistics.normalized_cross_moment = normalized_cross_moment;

    %% Pearson correlations
    minimum = min(J,[],2);
    maximum = max(J,[],2);
    constant_modes = minimum == maximum;

    fluctuation = J - first_moment;
    fluctuation(constant_modes,:) = 0;
    covariance = (fluctuation*fluctuation.')/N_membs;
    variance = diag(covariance);
    standard_deviation = sqrt(variance);
    pearson_correlations = covariance./(standard_deviation*standard_deviation.');

    statistics.covariance = covariance;
    statistics.variance = variance;
    statistics.standard_deviation = standard_deviation;
    statistics.pearson_correlations = pearson_correlations;

    %% Marginal fluctuation statistics
    third_central_moment = mean(fluctuation.^3,2);
    fourth_central_moment = mean(fluctuation.^4,2);

    statistics.fluctuation.mean_absolute = mean(abs(fluctuation),2);
    statistics.fluctuation.third_central_moment = third_central_moment;
    statistics.fluctuation.fourth_central_moment = fourth_central_moment;
    statistics.fluctuation.relative_variance = variance./first_moment.^2;
    statistics.fluctuation.coefficient_of_variation = standard_deviation./first_moment;
    statistics.fluctuation.fano_factor = variance./first_moment;
    statistics.fluctuation.g2 = second_moment./first_moment.^2;
    statistics.fluctuation.skewness = third_central_moment./standard_deviation.^3;
    statistics.fluctuation.kurtosis = fourth_central_moment./variance.^2;
    statistics.fluctuation.excess_kurtosis = statistics.fluctuation.kurtosis - 3;

    statistics.marginal.minimum = minimum;
    statistics.marginal.maximum = maximum;
    statistics.marginal.median = median(J,2);

    %% Marginal distributions
    N_bins = 100;
    marginal_edges = zeros(N_modes,N_bins + 1);
    marginal_centers = zeros(N_modes,N_bins);
    marginal_pdf = zeros(N_modes,N_bins);

    for mm = 1:N_modes
        [marginal_pdf(mm,:),marginal_edges(mm,:)] = ...
            histcounts(J(mm,:),N_bins,'Normalization','pdf');
        marginal_centers(mm,:) = ...
            (marginal_edges(mm,1:end-1) + marginal_edges(mm,2:end))/2;
    end

    statistics.marginal.N_bins = N_bins;
    statistics.marginal.edges = marginal_edges;
    statistics.marginal.centers = marginal_centers;
    statistics.marginal.fluctuation_edges = marginal_edges - first_moment;
    statistics.marginal.fluctuation_centers = marginal_centers - first_moment;
    statistics.marginal.pdf = marginal_pdf;

    %% Two-mode joint distributions
    N_joint_bins = min(50,max(10,floor(sqrt(N_membs/10))));
    N_pairs = N_modes*(N_modes - 1)/2;
    mode_pairs = zeros(N_pairs,2);
    joint_edges = zeros(N_modes,N_joint_bins + 1);
    joint_centers = zeros(N_modes,N_joint_bins);
    joint_pdf = zeros(N_joint_bins,N_joint_bins,N_pairs);

    for mm = 1:N_modes
        [~,joint_edges(mm,:)] = histcounts(J(mm,:),N_joint_bins);
        joint_centers(mm,:) = ...
            (joint_edges(mm,1:end-1) + joint_edges(mm,2:end))/2;
    end

    pp = 0;
    for mm = 1:N_modes
        for nn = mm+1:N_modes
            pp = pp + 1;
            mode_pairs(pp,:) = [mm,nn];
            joint_pdf(:,:,pp) = histcounts2( ...
                J(mm,:),J(nn,:),joint_edges(mm,:),joint_edges(nn,:), ...
                'Normalization','pdf');
        end
    end

    statistics.joint.N_bins = N_joint_bins;
    statistics.joint.N_pairs = N_pairs;
    statistics.joint.mode_pairs = mode_pairs;
    statistics.joint.edges = joint_edges;
    statistics.joint.centers = joint_centers;
    statistics.joint.pdf = joint_pdf;
end
