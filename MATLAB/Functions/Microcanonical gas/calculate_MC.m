function [T,mu,J] = calculate_MC(U,P,epsilon)
%CALCULATE_MC Microcanonical temperature, chemical potential, and means.
%   T and MU accept the same scalar or compatible array inputs as VOLUME.
%   The optional modal means J retain their scalar-state interface.

    epsilon = epsilon(:).';

    M = length(epsilon);

    %% Volume and its derivative
    V = volume(U,P,epsilon,[0 0]);
    dV = volume(U,P,epsilon,[1 0]);

    %% Temperature
    T = V./dV;

    %% Chemical potential
    mu = (U - (M - 2).*T)./P;

    %% First moment
    if nargout >= 3
        if ~isscalar(U) || ~isscalar(P)
            error('calculate_MC:ModalMeansBatch',...
                'The optional modal means require scalar U and P.')
        end

        J = zeros(M,1);
        for kk = 1:M
            alpha = zeros(1,M);
            alpha(kk) = 1;

            J(kk) = mixed_moments(U,P,epsilon,alpha);
        end
    end
end
