function [X,hierarchical] = system_hierarchy(X)
    hierarchical = isstruct(X) && isfield(X,'full');

    if ~hierarchical
        leaf = X;
        X = struct;
        X.full = leaf;
        X.sub = [];
    elseif ~isfield(X,'sub')
        X.sub = [];
    end
end
