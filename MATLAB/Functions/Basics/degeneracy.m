function [D,f,d] = degeneracy(epsilon)
    epsilon = epsilon(:).';

    tol = 1e-9;
    remaining = 1:length(epsilon);

    f = [];
    D = zeros(size(epsilon));
    d = cell(size(epsilon));

    while ~isempty(remaining)
        nn = remaining(1);
        d_n = remaining(abs(epsilon(remaining) - epsilon(nn)) <= tol);

        f(end+1) = nn; %#ok<AGROW>
        D(d_n) = length(d_n);

        for kk = d_n
            d{kk} = d_n;
        end

        remaining(ismember(remaining,d_n)) = [];
    end
end
