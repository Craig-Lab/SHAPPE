function stats = neighborhoodExposureStats(cells, types, r, minFocal)
% cells: Nx3 [x y type] with harmonized type labels
% types: 1xK desired order
% r: radius (e.g., 43)
% minFocal: minimum focal cells per row to keep (e.g., 15)

if nargin < 4, minFocal = 0; end
K = numel(types);

xy = cells(:,1:2);
tt = cells(:,3);

Mdl = KDTreeSearcher(xy);
nbrIdx = rangesearch(Mdl, xy, r);

typeToCol = containers.Map('KeyType','double','ValueType','double');
for k = 1:K
    typeToCol(types(k)) = k;
end

mu = zeros(K,K);      % mean counts
pi = zeros(K,K);      % mean fractions
e  = zeros(K,K);      % exposure P(N>=1)
nUsed = zeros(K,1);

% Loop through types
for ii = 1:K
    % Current type
    ti = types(ii);
    % Row numbers of cells of tpe ti
    idx_i = find(tt == ti);

    sumCounts = zeros(1,K);
    sumFrac   = zeros(1,K);
    sumExp    = zeros(1,K);
    used = 0;

    % Loop through cells of type ti
    for a = 1:numel(idx_i)
        % Row number of current cell
        p = idx_i(a);
        idx = nbrIdx{p};

        idx(idx==p) = [];

        % cell has no neighbor
        if(isempty(idx))
            counts = zeros(1,K);

            sumCounts = sumCounts + counts;
            sumFrac   = sumFrac   + counts;
            sumExp    = sumExp    + (counts > 0);

            used = used + 1;
            continue; 
        end

        % Types of neighbors
        tNbr = tt(idx);

        counts = zeros(1,K);
        % Loop through neighbors
        for q = 1:numel(tNbr)
            % Type of current neighbor
            tj = tNbr(q);
            if isKey(typeToCol, tj)
                jj = typeToCol(tj);
                counts(jj) = counts(jj) + 1;
            end
        end

        tot = sum(counts);

        sumCounts = sumCounts + counts;
        sumFrac   = sumFrac   + (counts / tot);
        sumExp    = sumExp    + (counts > 0);

        used = used + 1;
    end

    if used > 0
        mu(ii,:) = sumCounts / used;
        pi(ii,:) = sumFrac   / used;
        e(ii,:)  = sumExp    / used;
        nUsed(ii) = used;
    end
end

stats = struct('mu',mu,'pi',pi,'e',e,'nUsed',nUsed);
end