function [fvals,marks,weights] = ripleyK(coord, iType, jType, R, rvals,functionTransformation,shape)
    if(nargin<6)
        functionTransformation = "H";
        shape = "circle";
    elseif(nargin<7)
        shape = "circle";
    end

    rvals = rvals(:)';
    [rvals_sorted, ord] = sort(rvals);

    Xi = coord(coord(:,end)==iType, 1:2);
    Xj = coord(coord(:,end)==jType, 1:2);

    N_i = size(Xi,1);
    N_j = size(Xj,1);

    if N_i==0 || N_j==0
        fvals = nan(size(rvals));
        marks = [];
        return
    end

    if(strcmp(shape,"circle"))
        A = pi*R^2;
    elseif(strcmp(shape,"square"))
        A = (2.*R).^2;
    end
    lambda_i = N_i./A;
    if(iType==jType)
        lambda_j = (N_j-1)./A;
    else
        lambda_j = N_j./A;
    end

    norms_i = hypot(Xi(:,1), Xi(:,2));
    buffer = 0.001 * R;

    marks = [zeros(N_i,length(rvals_sorted))];

    % Create array of distances between points
    D = pdist2(Xi,Xj);

    K = zeros(1,length(rvals_sorted));

    weights(length(rvals_sorted)) = struct('wij',[]);
    % Loop through distances
    for rr=1:length(rvals_sorted)
        r = rvals_sorted(rr);
        wijMatrix = zeros(N_i,N_j);
        % Loop through points of type i
        for i=1:N_i
            m_i = 0;
            within_reach_point = find(D(i,:) <= r);
            % Remove itself if i=j
            if(iType==jType)
                within_reach_point = setdiff(within_reach_point,i);
            end
            if(~isempty(within_reach_point))
                % Make sure it's a row vector
                if(size(within_reach_point,1)>size(within_reach_point,2))
                    within_reach_point = within_reach_point';
                end
                xi = Xi(i,1);
                yi = Xi(i,2);
                for j=within_reach_point
                    dij = D(i,j);
                    if(strcmp(shape,"circle"))
                        li = norms_i(i);
                        boundary_condition = (li > 0) && (dij > 0) && (li + dij > (R + buffer));
                    elseif(strcmp(shape,"square"))
                        dist_to_edge = min([R-xi, xi+R, R-yi, yi+R]);
                        boundary_condition = dij > dist_to_edge;
                    end
                    winv = 1;
                    if(boundary_condition)
                        if(strcmp(shape,"circle"))
                            temp = (li^2 + dij^2 - R^2) / (2*li*dij);
                            temp = max(-1, min(1, temp));
                            phi  = 2*acos(temp);
                            winv = (2*pi) / phi;
                        elseif(strcmp(shape,"square"))
                            winv = isotropic_winv_square_mc(xi, yi, dij, R, 256);
                        end
                    end
                    wijMatrix(i,j) = winv;
                    m_i = m_i+winv;
                end
            end
            marks(i,rr) = (1./lambda_j).*m_i;
        end
        K(1,rr) = mean(marks(:,rr));
        weights(rr).wij = wijMatrix;
    end

    if(strcmp(functionTransformation,"L"))
        f_sorted = sqrt(K ./ pi);
    elseif(strcmp(functionTransformation,"H"))
        L = sqrt(K ./ pi);
        f_sorted = L - rvals_sorted;
    elseif(strcmp(functionTransformation,"K"))
        f_sorted = K;
    end
    fvals = zeros(size(rvals));
    fvals(ord) = f_sorted;

end