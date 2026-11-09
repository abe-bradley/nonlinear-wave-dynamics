function MM = mixed_moments(U,P,epsilon,alpha)
% MIXED_MOMENTS Normalized microcanonical moments of modal occupations.
%   The moment equation is accumulated in signed-log form and is retried
%   at variable precision only when double precision cannot resolve the
%   cancellation.
%
%     <prod J_k^alpha_k>
%         = prod(alpha_k!)*unnormalized moment/base volume.
%
%   The unnormalized moment uses the partial-fraction volume equation with
%   pole orders augmented by alpha. The original algorithm, which elucidates
%   the mathematical steps much clearer, is appended to the bottom of this
%   function as a comment. The more advanced algorithm used here only changes
%   the precision of evaluation.

    epsilon = epsilon(:).';
    alpha = alpha(:).';

    if ~isscalar(U) && ~isscalar(P) && ~isequal(size(U),size(P))
        error('mixed_moments:InputSize', ...
            'U and P must be scalars or arrays of equal size.')
    end

    [U_grid,P_grid,output_size] = expand_inputs(U,P);

    if all(alpha == 0)
        MM = ones(output_size);
        return
    end

    %% Degenerate levels and augmented pole orders
    [D,f,d] = degeneracy(epsilon);
    g = zeros(size(epsilon));

    for nn = f
        g(nn) = D(nn) + sum(alpha(d{nn}));
    end

    %% Exact spectral-edge limit
    [edge,edge_value] = spectral_edge_moment( ...
        U_grid,P_grid,epsilon,alpha,D,d);
    if all(edge,'all')
        MM = edge_value;
        return
    end

    %% Unnormalized numerator
    [sign_numerator,log_numerator] = unnormalized_moment_cached( ...
        U_grid,P_grid,epsilon,g,f,D);

    %% Normalize in the log domain
    g_base = zeros(size(epsilon));
    g_base(f) = D(f);
    [sign_V,log_V] = unnormalized_moment_cached( ...
        U_grid,P_grid,epsilon,g_base,f,D);
    log_prefactor = sum(gammaln(alpha + 1));

    sign_MM = sign_numerator.*sign_V;
    log_MM = log_prefactor + log_numerator - log_V;

    MM = zeros(output_size);
    nonzero = sign_MM ~= 0;
    MM(nonzero) = sign_MM(nonzero).*exp(log_MM(nonzero));
    MM(edge) = edge_value(edge);
end

function [U_grid,P_grid,output_size] = expand_inputs(U,P)
    if isscalar(U)
        output_size = size(P);
        U_grid = repmat(U,output_size);
    else
        output_size = size(U);
        U_grid = U;
    end

    if isscalar(P)
        P_grid = repmat(P,output_size);
    else
        P_grid = P;
    end
end

function [edge,value] = spectral_edge_moment( ...
        U,P,epsilon,alpha,D,d)
% At an exact spectral edge, all power lies in the extremal degenerate
% subspace. Its modal fractions are uniform on that simplex:
%
%   <prod J_k^alpha_k>
%       = P^A*Gamma(D)/Gamma(D+A)*prod Gamma(alpha_k+1).

    epsilon_min = min(epsilon);
    epsilon_max = max(epsilon);
    lower = epsilon_min.*P;
    upper = epsilon_max.*P;

    scale = max(abs(U),max(abs(lower),abs(upper)));
    tolerance = 32*eps(scale);

    near_lower = U < lower & U >= lower - tolerance;
    near_upper = U > upper & U <= upper + tolerance;
    U(near_lower) = lower(near_lower);
    U(near_upper) = upper(near_upper);

    lower_edge = isfinite(U) & isfinite(P) & P > 0 & U == lower;
    upper_edge = isfinite(U) & isfinite(P) & P > 0 ...
        & U == upper & ~lower_edge;

    edge = lower_edge | upper_edge;
    value = zeros(size(U));
    moment_order = sum(alpha);

    lower_index = find(epsilon == epsilon_min,1);
    lower_modes = d{lower_index};
    if all(alpha(setdiff(1:numel(epsilon),lower_modes)) == 0)
        log_coefficient = gammaln(D(lower_index)) ...
            - gammaln(D(lower_index) + moment_order) ...
            + sum(gammaln(alpha(lower_modes) + 1));
        value(lower_edge) = exp(log_coefficient ...
            + moment_order*log(P(lower_edge)));
    end

    upper_index = find(epsilon == epsilon_max,1);
    upper_modes = d{upper_index};
    if all(alpha(setdiff(1:numel(epsilon),upper_modes)) == 0)
        log_coefficient = gammaln(D(upper_index)) ...
            - gammaln(D(upper_index) + moment_order) ...
            + sum(gammaln(alpha(upper_modes) + 1));
        value(upper_edge) = exp(log_coefficient ...
            + moment_order*log(P(upper_edge)));
    end
end

function [sign_value,log_value] = unnormalized_moment_cached( ...
        U,P,epsilon,g,f,D)
% Cache the unnormalized integral by its augmented pole orders. Modal
% labels within one degenerate group produce the same g and therefore the
% same integral.

    persistent context g_cache sign_cache log_cache

    same_context = ~isempty(context) ...
        && isequaln(context.U,U) ...
        && isequaln(context.P,P) ...
        && isequaln(context.epsilon,epsilon);

    if ~same_context
        context = struct('U',U,'P',P,'epsilon',epsilon);
        g_cache = zeros(0,numel(g));
        sign_cache = {};
        log_cache = {};
    end

    cache_index = find(all(g_cache == g,2),1);
    if ~isempty(cache_index)
        sign_value = sign_cache{cache_index};
        log_value = log_cache{cache_index};

        % Keep frequently reused entries, especially the base volume, at
        % the back of the bounded cache.
        g_hit = g_cache(cache_index,:);
        g_cache(cache_index,:) = [];
        sign_cache(cache_index) = [];
        log_cache(cache_index) = [];
        g_cache(end+1,:) = g_hit;
        sign_cache{end+1} = sign_value;
        log_cache{end+1} = log_value;
        return
    end

    [sign_value,log_value] = unnormalized_moment( ...
        U,P,epsilon,g,f,D);

    g_cache(end+1,:) = g;
    sign_cache{end+1} = sign_value;
    log_cache{end+1} = log_value;

    cache_limit = 128;
    if size(g_cache,1) > cache_limit
        g_cache(1,:) = [];
        sign_cache(1) = [];
        log_cache(1) = [];
    end
end

function [sign_value,log_value] = unnormalized_moment( ...
        U,P,epsilon,g,f,D)

    output_size = size(U);
    [U,support] = physical_support(U,P,epsilon);

    if ~any(support,'all')
        sign_value = zeros(output_size);
        log_value = -Inf(output_size);
        return
    end

    %% Signed logarithms of the partial-fraction weights
    [~,sign_w,log_w,weights_reliable] = compute_weights(epsilon,g);

    %% Direct unnormalized moment equation
    log_positive = -Inf(output_size);
    log_negative = -Inf(output_size);
    g_tot = sum(g);

    for nn = f
        Qn = U - epsilon(nn).*P;
        Q = max(Qn,0);
        active = support & Qn >= 0;

        for mm = 1:g(nn)
            power_P = mm - 1;
            power_Q = g_tot - power_P - 2;

            coefficient_log = log_w(nn,mm) ...
                - gammaln(power_P + 1) - gammaln(power_Q + 1);

            [log_positive,log_negative] = accumulate_term( ...
                log_positive,log_negative,sign_w(nn,mm), ...
                coefficient_log,P,power_P,Q,power_Q,active);
        end
    end

    [~,sign_value,log_value,unresolved] = ...
        resolve_signed_logs(log_positive,log_negative);

    has_terms = isfinite(log_positive) | isfinite(log_negative);
    fallback = unresolved | (~weights_reliable & has_terms) ...
        | (has_terms & sign_value <= 0);

    if any(fallback,'all')
        g_base = zeros(size(g));
        g_base(f) = D(f);
        is_first_moment = isscalar(U) && isscalar(P) ...
            && sum(g - g_base) == 1 ...
            && all(D(f) == 1);

        if is_first_moment
            mode = find(g > g_base,1);
            [sign_all,log_all] = first_moments_vpa_cached(U,P,epsilon);
            sign_vpa = sign_all(mode);
            log_vpa = log_all(mode);
        else
            [sign_vpa,log_vpa] = unnormalized_moment_vpa( ...
                U(fallback),P(fallback),epsilon,g,f);
        end

        sign_value(fallback) = sign_vpa;
        log_value(fallback) = log_vpa;
    end
end

function [U,support] = physical_support(U,P,epsilon)
    lower = min(epsilon).*P;
    upper = max(epsilon).*P;

    scale = max(abs(U),max(abs(lower),abs(upper)));
    tolerance = 32*eps(scale);

    near_lower = U < lower & U >= lower - tolerance;
    near_upper = U > upper & U <= upper + tolerance;
    U(near_lower) = lower(near_lower);
    U(near_upper) = upper(near_upper);

    support = isfinite(U) & isfinite(P) & P > 0 ...
        & U >= lower & U <= upper;

    % Public normalized moments use the analytic edge limit above.
    support(U == lower | U == upper) = false;
end

function [log_positive,log_negative] = accumulate_term( ...
        log_positive,log_negative,coefficient_sign,coefficient_log, ...
        P,power_P,Q,power_Q,active)

    if coefficient_sign == 0 || coefficient_log == -Inf
        return
    end

    log_term = coefficient_log + zeros(size(P));
    term_sign = coefficient_sign*ones(size(P));
    valid = active;

    if power_P > 0
        valid = valid & P ~= 0;
        log_term(valid) = log_term(valid) ...
            + power_P*log(abs(P(valid)));
        term_sign(valid) = term_sign(valid).*sign(P(valid)).^power_P;
    end

    if power_Q > 0
        valid = valid & Q > 0;
        log_term(valid) = log_term(valid) ...
            + power_Q*log(Q(valid));
    end

    valid = valid & term_sign ~= 0 & ~isnan(log_term);
    positive = valid & term_sign > 0;
    negative = valid & term_sign < 0;

    log_positive(positive) = logaddexp( ...
        log_positive(positive),log_term(positive));
    log_negative(negative) = logaddexp( ...
        log_negative(negative),log_term(negative));
end

function [value,sign_value,log_value,unresolved] = ...
        resolve_signed_logs(log_positive,log_negative)

    value = zeros(size(log_positive));
    sign_value = zeros(size(log_positive));
    log_value = -Inf(size(log_positive));

    only_positive = isfinite(log_positive) & ~isfinite(log_negative);
    only_negative = ~isfinite(log_positive) & isfinite(log_negative);
    both = isfinite(log_positive) & isfinite(log_negative);

    sign_value(only_positive) = 1;
    log_value(only_positive) = log_positive(only_positive);

    sign_value(only_negative) = -1;
    log_value(only_negative) = log_negative(only_negative);

    positive_result = both & log_positive > log_negative;
    negative_result = both & log_negative > log_positive;

    sign_value(positive_result) = 1;
    log_value(positive_result) = logdiffexp( ...
        log_positive(positive_result),log_negative(positive_result));

    sign_value(negative_result) = -1;
    log_value(negative_result) = logdiffexp( ...
        log_negative(negative_result),log_positive(negative_result));

    nonzero = sign_value ~= 0;
    value(nonzero) = sign_value(nonzero).*exp(log_value(nonzero));

    separation = abs(log_positive - log_negative);
    unresolved = both & (separation < 1e-8 | ~isfinite(log_value));
    unresolved = unresolved | (nonzero & ~isfinite(value));
end

function value = logaddexp(log_x,log_y)
    value = log_x;

    use_y = log_x == -Inf;
    value(use_y) = log_y(use_y);

    y_positive_inf = isinf(log_y) & log_y > 0;
    value(y_positive_inf) = Inf;

    regular = isfinite(log_x) & isfinite(log_y);
    if any(regular,'all')
        largest = max(log_x(regular),log_y(regular));
        value(regular) = largest + log( ...
            exp(log_x(regular) - largest) ...
            + exp(log_y(regular) - largest));
    end
end

function value = logdiffexp(log_x,log_y)
    value = log_x + log(-expm1(log_y - log_x));
end

function [sign_value,log_value] = first_moments_vpa_cached(U,P,epsilon)
% All first moments share the simple-pole weights. Evaluating them together
% avoids repeating 105 nearly identical symbolic weight constructions.

    persistent cache

    if ~isempty(cache) ...
            && isequaln(cache.U,U) ...
            && isequaln(cache.P,P) ...
            && isequaln(cache.epsilon,epsilon)
        sign_value = cache.sign;
        log_value = cache.log;
        return
    end

    precision = max(64,ceil(numel(epsilon)/10));
    symbolic_value = converged_first_moments_vpa( ...
        U,P,epsilon,precision);

    sign_value = double(sign(symbolic_value));
    log_value = -Inf(size(symbolic_value));

    nonzero = sign_value ~= 0;
    if any(nonzero,'all')
        log_value(nonzero) = double(log(abs(symbolic_value(nonzero))));
    end

    cache = struct('U',U,'P',P,'epsilon',epsilon, ...
        'sign',sign_value,'log',log_value);
end

function result = converged_first_moments_vpa(U,P,epsilon,precision)
    if exist('vpa','file') ~= 2
        error('mixed_moments:VariablePrecisionUnavailable', ...
            'Adaptive cancellation recovery requires the Symbolic Math Toolbox.')
    end

    old_digits = digits;
    cleanup = onCleanup(@() digits(old_digits));

    previous = [];
    for attempt = 1:4
        digits(precision);
        current = first_moments_vpa_once(U,P,epsilon);

        if ~isempty(previous) && vpa_converged(current,previous)
            result = current;
            return
        end

        previous = current;
        precision = precision + 32;
    end

    error('mixed_moments:Precision', ...
        'The first moments did not converge at variable precision.')
end

function unnormalized_moments = first_moments_vpa_once(U,P,epsilon)
% For alpha=e_k and distinct energies, all augmented weights follow from
% the simple-pole weights. This is the same partial-fraction equation,
% evaluated for every k in one matrix operation.

    epsilon_vpa = vpa(sym(epsilon(:),'f'));
    M = numel(epsilon);

    base_weights = compute_weights(epsilon_vpa,ones(M,1));
    base_weights = base_weights(:,1);

    Qn = U - epsilon(:)*P;
    U_vpa = vpa(sym(U,'f'));
    P_vpa = vpa(sym(P,'f'));
    Q = (U_vpa - epsilon_vpa*P_vpa).*vpa(Qn > 0);

    A = Q.^(M-1)/gamma(vpa(M));
    B = P_vpa*Q.^(M-2)/gamma(vpa(M-1));

    energy_difference = epsilon_vpa - epsilon_vpa.';
    energy_difference(1:M+1:end) = vpa(1);

    inverse_difference = 1./energy_difference;
    inverse_difference(1:M+1:end) = vpa(0);

    s1 = sum(inverse_difference,2);
    cross = inverse_difference*(base_weights.*A);

    unnormalized_moments = ...
        cross + base_weights.*s1.*A + base_weights.*B;
end

function [sign_value,log_value] = unnormalized_moment_vpa( ...
        U,P,epsilon,g,f)

    precision = max(64,ceil(sum(g)/10));
    symbolic_value = converged_unnormalized_moment_vpa( ...
        U,P,epsilon,g,f,precision);

    sign_value = double(sign(symbolic_value));
    log_value = -Inf(size(symbolic_value));

    nonzero = sign_value ~= 0;
    if any(nonzero,'all')
        log_value(nonzero) = double(log(abs(symbolic_value(nonzero))));
    end
end

function result = converged_unnormalized_moment_vpa( ...
        U,P,epsilon,g,f,precision)

    if exist('vpa','file') ~= 2
        error('mixed_moments:VariablePrecisionUnavailable', ...
            'Adaptive cancellation recovery requires the Symbolic Math Toolbox.')
    end

    old_digits = digits;
    cleanup = onCleanup(@() digits(old_digits));

    previous = [];
    for attempt = 1:4
        digits(precision);
        current = unnormalized_moment_vpa_once(U,P,epsilon,g,f);

        if ~isempty(previous) && vpa_converged(current,previous)
            result = current;
            return
        end

        previous = current;
        precision = precision + 32;
    end

    error('mixed_moments:Precision', ...
        'The mixed moment did not converge at variable precision.')
end

function result = unnormalized_moment_vpa_once(U,P,epsilon,g,f)
    U_vpa = vpa(sym(U,'f'));
    P_vpa = vpa(sym(P,'f'));
    epsilon_vpa = vpa(sym(epsilon,'f'));
    w = compute_weights(epsilon_vpa,g);

    result = sym(zeros(size(U_vpa)));
    g_tot = sum(g);

    for nn = f
        Qn = U - epsilon(nn).*P;
        active = Qn >= 0;
        Q = (U_vpa - epsilon_vpa(nn).*P_vpa).*vpa(Qn > 0);

        for mm = 1:g(nn)
            power_P = mm - 1;
            power_Q = g_tot - power_P - 2;

            result = result + w(nn,mm) ...
                .* P_vpa.^power_P/gamma(vpa(power_P + 1)) ...
                .* Q.^power_Q/gamma(vpa(power_Q + 1)) ...
                .* vpa(active);
        end
    end

end

function converged = vpa_converged(current,previous)
    converged = true;

    for kk = 1:numel(current)
        if isAlways(current(kk) == 0)
            converged = converged && isAlways(previous(kk) == 0);
        else
            relative_change = abs((current(kk) - previous(kk))/current(kk));
            converged = converged && double(relative_change) < 1e-12;
        end
    end
end



%% ARCHIVED MIXED MOMENTS FUNCTION
% function MM = mixed_moments(U,P,epsilon,alpha)
%     epsilon = epsilon(:).';
%     alpha = alpha(:).'; % non-negative integer

%     %% Sift first occurrences
%     [D,f,d] = degeneracy(epsilon);
%     g = zeros(size(epsilon));

%     for nn = f
%         g(nn) = D(nn) + sum(alpha(d{nn}));
%     end

%     g_tot = sum(g);

%     %% Weights: w_{m,n}[g]
%     w = compute_weights(epsilon,g);

%     %% Main sum
%     Phi = 0;
%     for nn = f
%         Q = max(U - epsilon(nn)*P, 0);

%         for mm = 1:g(nn)
%             ii = mm - 1;
%             jj = g_tot - ii - 2;

%             Phi = Phi + w(nn,mm)*P^ii/factorial(ii) *...
%                 Q^jj/factorial(jj);
%         end
%     end

%     %% Normalize
%     pre = prod(factorial(alpha));
%     V = volume(U,P,epsilon,[0 0]);

%     MM = pre*Phi/V;
% end