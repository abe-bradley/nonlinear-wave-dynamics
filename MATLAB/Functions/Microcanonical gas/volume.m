function [dV,sign_dV,log_dV] = volume(U,P,epsilon,dO)
% VOLUME Derivatives of the constrained microcanonical phase-space volume.
%   The defining partial-fraction sum is evaluated in signed-log form.
%   Entries with unresolved cancellation are recomputed adaptively using
%   variable precision. Entries outside the physical energy interval and
%   known-zero spectral boundaries are returned exactly without summation.
%
%   Before evaluation, the energy origin is shifted to the spectral mean,
%
%     epsilon -> epsilon - mean(epsilon),
%     U       -> U - mean(epsilon)*P.
%
%   This leaves every constraint U - epsilon_n*P invariant. Centered
%   positive energies are then evaluated through the exact reflection
%
%     V(U,P;epsilon) = V(-U,P;-epsilon).
%
%   These analytic coordinate symmetries require no symmetry of the
%   underlying spectrum. Derivative factors retain the original energy
%   coordinates, so outputs remain derivatives with respect to the input
%   U and P.
%
%   The main loop is the direct equation
%
%     sum w_nm*C*(-epsilon_n)^p
%         * P^(i-I)/(i-I)! * [U-epsilon_n P]_+^(j-J)/(j-J)! .
%
%   The original algorithm, which elucidates the mathematical steps much 
%   clearer, is appended to the bottom of this function as a comment. The
%   more advanced algorithm used here only changes the precision of evaluation.

    epsilon_original = epsilon(:).';
    dO = dO(:).';

    if ~isscalar(U) && ~isscalar(P) && ~isequal(size(U),size(P))
        error('volume:InputSize',...
            'U and P must be scalars or arrays of equal size.')
    end

    [U_grid,P_grid,output_size] = expand_inputs(U,P);

    %% Center and reflect the energy coordinates
    energy_center = mean(epsilon_original);
    epsilon_centered = epsilon_original - energy_center;
    U_centered = U_grid - energy_center.*P_grid;
    reflect = U_centered > 0;

    dV = zeros(output_size);
    sign_dV = zeros(output_size);
    log_dV = -Inf(output_size);

    direct = ~reflect;
    if any(direct,'all')
        [dV(direct),sign_dV(direct),log_dV(direct)] = ...
            volume_centered(U_centered(direct),P_grid(direct),...
            epsilon_centered,epsilon_original,dO);
    end

    if any(reflect,'all')
        [dV(reflect),sign_dV(reflect),log_dV(reflect)] = ...
            volume_centered(-U_centered(reflect),P_grid(reflect),...
            -epsilon_centered,-epsilon_original,dO);

        if mod(dO(1),2) == 1
            dV(reflect) = -dV(reflect);
            sign_dV(reflect) = -sign_dV(reflect);
        end
    end
end

function [dV,sign_dV,log_dV] = volume_centered( ...
        U_grid,P_grid,epsilon,epsilon_derivative,dO)
    output_size = size(U_grid);

    %% Degenerate energy levels
    [D,f,~] = degeneracy(epsilon);
    g = zeros(size(epsilon));
    g(f) = D(f);
    g_tot = sum(g);

    %% Physical support
    [U_grid,support] = physical_support( ...
        U_grid,P_grid,epsilon,D,dO);

    if ~any(support,'all')
        dV = zeros(output_size);
        sign_dV = zeros(output_size);
        log_dV = -Inf(output_size);
        return
    end

    %% Signed logarithms of the partial-fraction weights
    [~,sign_w,log_w,weights_reliable] = ...
        compute_weights(epsilon,g,true);

    %% Defining volume equation
    log_positive = -Inf(output_size);
    log_negative = -Inf(output_size);

    dO_U = dO(1);
    dO_P = dO(2);

    for nn = f
        Qn = U_grid - epsilon(nn).*P_grid;
        Q = max(Qn,0);

        for mm = 1:g(nn)
            ii = mm - 1;
            jj = g_tot - ii - 2;

            for pp = 0:dO_P
                II = dO_P - pp;
                JJ = dO_U + pp;

                if II > ii || JJ > jj
                    continue
                end

                power_P = ii - II;
                power_Q = jj - JJ;

                [factor_sign,factor_log] = derivative_factor( ...
                    dO_P,pp,epsilon_derivative(nn));
                coefficient_sign = sign_w(nn,mm)*factor_sign;
                coefficient_log = log_w(nn,mm) + factor_log ...
                    - gammaln(power_P + 1) - gammaln(power_Q + 1);

                if JJ == 0
                    active = support & Qn >= 0;
                else
                    active = support & Qn > 0;
                end

                [log_positive,log_negative] = accumulate_term( ...
                    log_positive,log_negative,coefficient_sign, ...
                    coefficient_log,P_grid,power_P,Q,power_Q,active);
            end
        end
    end

    [dV,sign_dV,log_dV,unresolved] = ...
        resolve_signed_logs(log_positive,log_negative);

    has_terms = isfinite(log_positive) | isfinite(log_negative);
    fallback = unresolved | (~weights_reliable & has_terms);

    % The undifferentiated volume is nonnegative. Other derivatives are
    % allowed to carry either sign.
    if isequal(dO,[0 0])
        fallback = fallback | (has_terms & sign_dV <= 0);
    end

    if any(fallback,'all')
        [dV_vpa,sign_vpa,log_vpa] = volume_vpa_cached( ...
            U_grid(fallback),P_grid(fallback),epsilon,...
            epsilon_derivative,dO,g,f);

        dV(fallback) = dV_vpa;
        sign_dV(fallback) = sign_vpa;
        log_dV(fallback) = log_vpa;
    end
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

function [U,support] = physical_support(U,P,epsilon,D,dO)
% The constrained simplex is empty unless
%
%   epsilon_min*P <= U <= epsilon_max*P,   P > 0.
%
% Snap roundoff-sized exterior excursions to the exact boundary. At a
% spectral edge the volume vanishes with power M-D_edge-1, so derivatives
% below that order are known zeros and need no numerical cancellation.

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

    support = isfinite(U) & isfinite(P) & P > 0 ...
        & U >= lower & U <= upper;

    derivative_order = sum(dO);
    lower_index = find(epsilon == epsilon_min,1);
    upper_index = find(epsilon == epsilon_max,1);
    lower_power = numel(epsilon) - D(lower_index) - 1;
    upper_power = numel(epsilon) - D(upper_index) - 1;

    known_zero = U == lower & lower_power > derivative_order;
    known_zero = known_zero ...
        | (U == upper & upper_power > derivative_order);
    support(known_zero) = false;
end

function [factor_sign,factor_log] = derivative_factor(dO_P,pp,epsilon_n)
    factor_log = gammaln(dO_P + 1) ...
        - gammaln(pp + 1) - gammaln(dO_P - pp + 1);

    if pp == 0
        factor_sign = 1;
    elseif epsilon_n == 0
        factor_sign = 0;
        factor_log = -Inf;
    else
        factor_sign = sign(-epsilon_n)^pp;
        factor_log = factor_log + pp*log(abs(epsilon_n));
    end
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

function [value,sign_value,log_value] = volume_vpa_cached( ...
        U,P,epsilon,epsilon_derivative,dO,g,f)
% Cache the most recent expensive fallback; calculate_MC requests the same
% base volume once for every modal moment.

    persistent cache

    for kk = 1:numel(cache)
        if isequaln(cache(kk).U,U) ...
                && isequaln(cache(kk).P,P) ...
                && isequaln(cache(kk).epsilon,epsilon) ...
                && isequaln(cache(kk).epsilon_derivative,epsilon_derivative) ...
                && isequaln(cache(kk).dO,dO)
            value = cache(kk).value;
            sign_value = cache(kk).sign;
            log_value = cache(kk).log;
            return
        end
    end

    precision = max(64,ceil(sum(g)/10));
    symbolic_value = converged_volume_vpa( ...
        U,P,epsilon,epsilon_derivative,dO,g,f,precision);

    sign_value = double(sign(symbolic_value));
    log_value = -Inf(size(symbolic_value));

    nonzero = sign_value ~= 0;
    if any(nonzero,'all')
        log_value(nonzero) = double(log(abs(symbolic_value(nonzero))));
    end
    value = double(symbolic_value);

    entry = struct('U',U,'P',P,'epsilon',epsilon,...
        'epsilon_derivative',epsilon_derivative,'dO',dO, ...
        'value',value,'sign',sign_value,'log',log_value);
    cache = [entry,cache];
    cache = cache(1:min(2,numel(cache)));
end

function result = converged_volume_vpa( ...
        U,P,epsilon,epsilon_derivative,dO,g,f,precision)
    if exist('vpa','file') ~= 2
        error('volume:VariablePrecisionUnavailable', ...
            'Adaptive cancellation recovery requires the Symbolic Math Toolbox.')
    end

    old_digits = digits;
    cleanup = onCleanup(@() digits(old_digits));

    previous = [];
    for attempt = 1:4
        digits(precision);
        current = volume_vpa_once( ...
            U,P,epsilon,epsilon_derivative,dO,g,f);

        if ~isempty(previous) && vpa_converged(current,previous)
            result = current;
            return
        end

        previous = current;
        precision = precision + 32;
    end

    error('volume:Precision', ...
        'The phase-space volume did not converge at variable precision.')
end

function result = volume_vpa_once( ...
        U,P,epsilon,epsilon_derivative,dO,g,f)
    U_vpa = vpa(sym(U,'f'));
    P_vpa = vpa(sym(P,'f'));
    epsilon_vpa = vpa(sym(epsilon,'f'));
    epsilon_derivative_vpa = vpa(sym(epsilon_derivative,'f'));
    w = compute_weights(epsilon_vpa,g);

    result = sym(zeros(size(U_vpa)));
    g_tot = sum(g);
    dO_U = dO(1);
    dO_P = dO(2);

    for nn = f
        Qn = U - epsilon(nn).*P;
        Q = (U_vpa - epsilon_vpa(nn).*P_vpa).*vpa(Qn > 0);

        for mm = 1:g(nn)
            ii = mm - 1;
            jj = g_tot - ii - 2;

            for pp = 0:dO_P
                II = dO_P - pp;
                JJ = dO_U + pp;

                if II > ii || JJ > jj
                    continue
                end

                power_P = ii - II;
                power_Q = jj - JJ;

                if JJ == 0
                    active = Qn >= 0;
                else
                    active = Qn > 0;
                end

                term = w(nn,mm)*nchoosek(dO_P,pp) ...
                    * (-epsilon_derivative_vpa(nn))^pp ...
                    .* P_vpa.^power_P/gamma(vpa(power_P + 1)) ...
                    .* Q.^power_Q/gamma(vpa(power_Q + 1)) ...
                    .* vpa(active);

                result = result + term;
            end
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


%% ARCHIVED VOLUME FUNCTION
% function dV = volume(U,P,epsilon,dO)
%     epsilon = epsilon(:).';
%     dO = dO(:).';

%     if ~isscalar(U) && ~isscalar(P) && ~isequal(size(U),size(P))
%         error('volume:InputSize',...
%             'U and P must be scalars or arrays of equal size.')
%     end

%     if isscalar(U)
%         output_size = size(P);
%     else
%         output_size = size(U);
%     end

%     %% Sift first occurrences
%     [D,f,~] = degeneracy(epsilon);
%     g = zeros(size(epsilon));
%     g(f) = D(f);

%     g_tot = sum(g);

%     %% Weights: w_{m,n}[g]
%     w = compute_weights(epsilon,g);

%     %% Main sum
%     dV = zeros(output_size);
%     dO_U = dO(1);
%     dO_P = dO(2);

%     for nn = f                              % First representative level index n
%         Qn = U - epsilon(nn).*P;            % Positive-part argument Q_n = U - epsilon_n P

%         for mm = 1:g(nn)                    % Pole-order summation index m
%             ii = mm - 1;                    % Power of P in the rewritten summand
%             jj = g_tot - ii - 2;            % Power of Q_n in the rewritten summand

%             for pp = 0:dO_P               % Number of P derivatives applied to Q_n
%                 II = dO_P - pp;           % Number of P derivatives applied to P^ii
%                 JJ = dO_U + pp;           % Total number of derivatives applied to Q_n^jj

%                 if II <= ii && JJ <= jj     % Skip terms killed by too many derivatives
%                     if JJ == 0              % No derivative hits the positive-part factor
%                         Q = max(Qn, 0);     % Use Q_n^+ for the undifferentiated factor
%                     else                    % Derivative hits the positive-part factor
%                         active = Qn > 0;     % Strictly positive entries contribute
%                         Q = max(Qn, 0);     % Inactive entries are masked below
%                     end

%                     term = w(nn,mm)*nchoosek(dO_P,pp)*(-epsilon(nn))^pp .*...
%                         P.^(ii - II)./factorial(ii - II) .*...
%                         Q.^(jj - JJ)./factorial(jj - JJ);

%                     if JJ > 0
%                         term(~active) = 0;  % Avoid inactive 0^0 contributions
%                     end

%                     dV = dV + term;
%                 end
%             end
%         end
%     end
% end