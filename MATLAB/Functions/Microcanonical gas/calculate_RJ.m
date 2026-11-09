function [T,mu,J] = calculate_RJ(U,P,epsilon)
    % Uses the convention that U = epsilon*J0.
    % If the energy is given as U = -beta*J0, then epsilon = -beta must be used.

    [epsilon,I_order] = sort(epsilon,'descend');

    M = length(epsilon);
    rho = U - P*epsilon;

    interval_side = sign(sum(rho));
    delta = 200;
    if interval_side == -1
        T0 = rho(M) + delta;
    else
        T0 = rho(1) - delta;
    end

    T = fzero(@f,T0);
    mu = U/P - M*T/P;

    if nargout >= 3
        J = T./(epsilon - mu);
        J(I_order) = J;
    end
    
    function y = f(T)
        y = sum(T./(M*T - rho)) - 1;

        if interval_side == -1
            y = y.*(sign(M*T - rho(M)) + 1)/2 + abs((sign(M*T - rho(M)) - 1)/2);
        else
            y = y.*(-sign(M*T - rho(1)) + 1)/2 + abs((-sign(M*T - rho(1)) - 1)/2);
        end
    end
end
