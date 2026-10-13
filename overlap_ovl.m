function ovl = overlap_ovl(x1, x2, nGrid)
    if nargin < 3, nGrid = 1024; end
    x1 = x1(:); x2 = x2(:);

    lo = min([x1; x2]);
    hi = max([x1; x2]);
    xi = linspace(lo, hi, nGrid);

    f1 = ksdensity(x1, xi, 'Function','pdf');
    f2 = ksdensity(x2, xi, 'Function','pdf');

    % Numerical integration via the trapezoidal method
    ovl = trapz(xi, min(f1,f2));
    ovl = max(0, min(1, ovl));
end