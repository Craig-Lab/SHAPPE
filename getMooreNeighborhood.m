function neigh_bins_idx = getMooreNeighborhood(nBin,b,layer)
    % Center bin (row/col)
    [i0, j0] = ind2sub([nBin, nBin], b);

    % Ring offsets (Chebyshev distance = layer)
    [di, dj] = ndgrid(-layer:layer, -layer:layer);
    ringMask = max(abs(di), abs(dj)) == layer;

    rr = i0 + di(ringMask);
    cc = j0 + dj(ringMask);

    % --- 1) keep only in-bounds first ---
    keepBounds = (rr >= 1 & rr <= nBin & cc >= 1 & cc <= nBin);
    rr = rr(keepBounds);
    cc = cc(keepBounds);

    if isempty(rr)
        neigh_bins_idx = [];
        return
    end

    % --- 2) build inside-circle mask on the grid ---
    % c = (nBin+1)/2;
    % [xi, yi] = ndgrid(1:nBin, 1:nBin);
    % X = xi - c;
    % Y = yi - c;
    % insideCircle = (X.^2 + Y.^2 <= (nBin/2)^2);

    % --- 3) now safe to linear-index and filter by circle ---
    % lin = sub2ind([nBin, nBin], rr, cc);
    % keepCircle = insideCircle(lin);
    % 
    % rr = rr(keepCircle);
    % cc = cc(keepCircle);

    % Final neighbor bin indices
    neigh_bins_idx = sub2ind([nBin, nBin], rr, cc);
    neigh_bins_idx = sort(neigh_bins_idx(:))';
end