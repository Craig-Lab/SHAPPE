function winv = isotropic_winv_square_mc(x, y, rvec, R, nAng)
% Isotropic edge correction for square [-R,R]^2 via angle sampling.
% Deterministic angles => reproducible and less noisy.

    rvec = rvec(:);
    winv = ones(size(rvec));

    theta = linspace(0,2*pi,nAng+1)'; 
    theta(end) = [];
    ct = cos(theta);
    st = sin(theta);

    for k = 1:numel(rvec)
        r = rvec(k);
        if r <= 0
            winv(k) = 1;
            continue;
        end

        X = x + r*ct;
        Y = y + r*st;

        inside = (abs(X) <= R) & (abs(Y) <= R);
        p = mean(inside);
        p = max(p, 1/nAng);   % safety

        winv(k) = 1/p;
    end
end
