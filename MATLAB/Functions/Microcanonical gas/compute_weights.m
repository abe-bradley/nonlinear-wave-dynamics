function [w,sign_w,log_w,reliable] = compute_weights( ...
        epsilon,g,recover_unreliable)
%COMPUTE_WEIGHTS Partial-fraction weights for degenerate energy levels.
%   W = COMPUTE_WEIGHTS(EPSILON,G) preserves the original numeric API.
%   Internal callers may also request SIGN_W and LOG_W so signed sums can
%   be accumulated without first expanding very large weight magnitudes.
%   The repeated-pole coefficient recurrence is unchanged; its additions
%   are simply performed as signed logarithms. Internal volume evaluations
%   may request adaptive variable-precision recovery of an unreliable
%   weight matrix before returning to double-precision signed-log
%   arithmetic.

    if nargin < 3
        recover_unreliable = false;
    end

    epsilon = epsilon(:);
    g = g(:);

    if isa(epsilon,'sym')
        w = compute_weights_vpa(epsilon,g);
        sign_w = [];
        log_w = [];
        reliable = true;
        return
    end

    %% Unique energy levels
    [~,f,~] = degeneracy(epsilon);
    epsilon_f = epsilon(f);
    g_f = g(f);

    M = length(f);
    g_max = max(g_f);

    %% Pairwise energy differences
    Delta = epsilon_f - epsilon_f.';
    Delta(1:M+1:end) = 1;

    sign_Delta = sign(Delta);
    log_Delta = log(abs(Delta));

    %% Logarithmic power sums
    sign_s = zeros(M,g_max-1);
    log_s = -Inf(M,g_max-1);
    log_g = log(g_f);
    reliable = true;

    for pp = 1:g_max-1
        sign_terms = ((-1)^pp)*sign_Delta.^pp;
        log_terms = log_g - pp*log_Delta;

        sign_terms(1:M+1:end) = 0;
        log_terms(1:M+1:end) = -Inf;

        for nn = 1:M
            [sign_s(nn,pp),log_s(nn,pp),resolved] = ...
                signed_logsumexp(sign_terms(:,nn),log_terms(:,nn));
            reliable = reliable && resolved;
        end
    end

    %% Leading coefficient
    L = log(abs(Delta));
    L(1:M+1:end) = 0;
    log_c0 = -(L.'*g_f);

    negative = epsilon_f < epsilon_f.';
    negative(1:M+1:end) = false;
    parity = mod(negative.'*g_f,2);
    sign_c0 = 1 - 2*parity;

    %% Repeated-pole coefficients
    sign_c = zeros(M,g_max);
    log_c = -Inf(M,g_max);
    sign_c(:,1) = sign_c0;
    log_c(:,1) = log_c0;

    for kk = 0:g_max-2
        for nn = 1:M
            sign_terms = sign_c(nn,kk+1:-1:1).*sign_s(nn,1:kk+1);
            log_terms = log_c(nn,kk+1:-1:1) + log_s(nn,1:kk+1);

            [sign_c(nn,kk+2),log_c(nn,kk+2),resolved] = ...
                signed_logsumexp(sign_terms,log_terms);
            reliable = reliable && resolved;
            log_c(nn,kk+2) = log_c(nn,kk+2) - log(kk + 1);
        end
    end

    %% Restore the original weight layout
    sign_w = zeros(numel(epsilon),g_max);
    log_w = -Inf(numel(epsilon),g_max);

    for nn = 1:M
        gm = g_f(nn);
        sign_w(f(nn),1:gm) = sign_c(nn,gm:-1:1);
        log_w(f(nn),1:gm) = log_c(nn,gm:-1:1);
    end

    if ~reliable && recover_unreliable
        [sign_vpa,log_vpa,recovered] = ...
            converged_weights_vpa(epsilon,g);

        if recovered
            sign_w = sign_vpa;
            log_w = log_vpa;
            reliable = true;
        end
    end

    w = zeros(size(sign_w));
    active = sign_w ~= 0;
    w(active) = sign_w(active).*exp(log_w(active));
end

function [sign_w,log_w,recovered] = converged_weights_vpa(epsilon,g)
    sign_w = [];
    log_w = [];
    recovered = false;

    if exist('vpa','file') ~= 2
        return
    end

    old_digits = digits;
    cleanup = onCleanup(@() digits(old_digits));

    precision = max(64,ceil(sum(g)/10));
    previous = [];

    active = false(numel(epsilon),max(g));
    f = find(g > 0);
    for nn = f.'
        active(nn,1:g(nn)) = true;
    end
    active_indices = find(active);

    for attempt = 1:4
        digits(precision);
        epsilon_vpa = vpa(sym(epsilon,'f'));
        current = compute_weights_vpa(epsilon_vpa,g);

        if ~isempty(previous) && weights_vpa_converged( ...
                current,previous,active_indices)
            sign_w = double(sign(current));
            log_w = -Inf(size(sign_w));

            nonzero = sign_w ~= 0;
            if any(nonzero,'all')
                log_w(nonzero) = double(log(abs(current(nonzero))));
            end
            recovered = true;
            return
        end

        previous = current;
        precision = precision + 32;
    end
end

function w = compute_weights_vpa(epsilon,g)
% Direct variable-precision form used only by an adaptive fallback.

    f = find(g > 0);
    epsilon_f = epsilon(f);
    g_f = g(f);

    M = length(f);
    g_max = max(g_f);

    Delta = epsilon_f - epsilon_f.';
    Delta(1:M+1:end) = vpa(1);

    inv_Delta = 1./Delta;
    inv_Delta(1:M+1:end) = vpa(0);

    s = sym(zeros(M,g_max-1));
    inv_Delta_p = inv_Delta;
    for pp = 1:g_max-1
        s(:,pp) = ((-1)^pp)*(inv_Delta_p.'*vpa(g_f));
        inv_Delta_p = inv_Delta_p.*inv_Delta;
    end

    c = sym(zeros(M,g_max));
    for nn = 1:M
        c(nn,1) = vpa(1);
        for ll = 1:M
            if ll ~= nn
                c(nn,1) = c(nn,1)*Delta(ll,nn)^(-g_f(ll));
            end
        end
    end

    for kk = 0:g_max-2
        c(:,kk+2) = sum(c(:,kk+1:-1:1).*s(:,1:kk+1),2)/(kk + 1);
    end

    w = sym(zeros(numel(epsilon),g_max));
    for nn = 1:M
        gm = g_f(nn);
        w(f(nn),1:gm) = c(nn,gm:-1:1);
    end
end

function converged = weights_vpa_converged( ...
        current,previous,active_indices)
    converged = true;

    for kk = active_indices.'
        current_zero = isAlways(current(kk) == 0);
        previous_zero = isAlways(previous(kk) == 0);

        if current_zero || previous_zero
            converged = current_zero && previous_zero;
        else
            relative_change = abs( ...
                (current(kk) - previous(kk))/current(kk));
            converged = double(relative_change) < 1e-12;
        end

        if ~converged
            return
        end
    end
end

function [sign_sum,log_sum,reliable] = signed_logsumexp(sign_terms,log_terms)
    log_positive = logsumexp(log_terms(sign_terms > 0));
    log_negative = logsumexp(log_terms(sign_terms < 0));
    reliable = true;

    if isfinite(log_positive) && isfinite(log_negative)
        reliable = abs(log_positive - log_negative) >= 1e-8;
    end

    if log_positive == -Inf && log_negative == -Inf
        sign_sum = 0;
        log_sum = -Inf;
    elseif log_positive == log_negative
        sign_sum = 0;
        log_sum = -Inf;
    elseif log_positive > log_negative
        sign_sum = 1;
        log_sum = logdiffexp(log_positive,log_negative);
    else
        sign_sum = -1;
        log_sum = logdiffexp(log_negative,log_positive);
    end
end

function value = logsumexp(log_terms)
    log_terms = log_terms(isfinite(log_terms));

    if isempty(log_terms)
        value = -Inf;
    else
        largest = max(log_terms);
        value = largest + log(sum(exp(log_terms - largest)));
    end
end

function value = logdiffexp(log_x,log_y)
    if log_y == -Inf
        value = log_x;
    elseif log_x == log_y
        value = -Inf;
    else
        value = log_x + log1p(-exp(log_y - log_x));
    end
end
