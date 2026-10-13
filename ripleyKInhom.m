function [fvals,marks] = ripleyKInhom(coord, iType, jType, R, rvals,h,functionTransformation,shape)
    if(nargin<7)
        functionTransformation = "H";
        shape = "circle";
    elseif(nargin<8)
        shape = "circle";
    end

    rvals = rvals(:)';
    [rvals_sorted, ord] = sort(rvals);

    % Positions of cells of type i
    Xi = coord(coord(:,end)==iType, 1:2);
    % Positions of cells of type j
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
    % Generate a polyshape for the window
    if(strcmp(shape,"circle"))
        theta = linspace(0,2.*pi,300);
        x = R.*cos(theta);
        y = R.*sin(theta);
        ROI = polyshape(x,y);
    elseif(strcmp(shape,"square"))
        ROI = polyshape([-R -R R R],[R -R -R R]);
    end

    norms_i = hypot(Xi(:,1), Xi(:,2));
    buffer = 0.001 * R;

    marks = [zeros(N_i,length(rvals_sorted))];

    % Create array of distances between points
    D = pdist2(Xi,Xj);

    % Estimate local densities
    %%% Type i
    Di = squareform(pdist(Xi));
    lambda_i = zeros(N_i,1);
    for i=1:N_i
        pos = Xi(i,1:2);
        AreaMovingWindow = pi.*h.^2;
        normaliseArea = false;
        if(strcmp(shape,"square"))
            if(abs(abs(pos(1,1))-R)<h || abs(abs(pos(1,2))-R)<h)
                normaliseArea = true;
            end
        elseif(strcmp(shape,"circle"))
            if(pos(1,1).^2+pos(1,2).^2>(R-h).^2)
                normaliseArea = true;
            end
        end
        if(normaliseArea)
            theta = linspace(0,2.*pi,300);
            x_movingWindow = pos(1,1)+h.*cos(theta);
            y_movingWindow = pos(1,2)+h.*sin(theta);
            movingWindow = polyshape(x_movingWindow,y_movingWindow);
            AreaMovingWindow = area(intersect(ROI,movingWindow));
        end
        distances_to_other_points = Di(i,:);
        % Removes points to far
        distances_to_other_points(i) = [];
        distances_to_other_points(distances_to_other_points>h) = [];
        % Reweight points within reach
        distances_to_other_points = (3./(4.*h)).*(1-(distances_to_other_points.^2)./(h.^2));
        % lambda_i(i,1) = sum(distances_to_other_points)./AreaMovingWindow;
        lambda_i(i,1) = numel(distances_to_other_points)./AreaMovingWindow;
    end
    if(isempty(lambda_i(lambda_i(:,1)>0,1)))
        lambda_i(:) = N_i./A;
    else
        lambda_i(lambda_i(:,1)==0,1) = min(lambda_i(lambda_i(:,1)>0,1));
    end
    %%% Type j
    Dj = squareform(pdist(Xj));
    lambda_j = zeros(N_j,1);
    for i=1:N_j
        pos = Xj(i,1:2);
        AreaMovingWindow = pi.*h.^2;
        normaliseArea = false;
        if(strcmp(shape,"square"))
            if(abs(abs(pos(1,1))-R)<h || abs(abs(pos(1,2))-R)<h)
                normaliseArea = true;
            end
        elseif(strcmp(shape,"circle"))
            if(pos(1,1).^2+pos(1,2).^2>(R-h).^2)
                normaliseArea = true;
            end
        end
        if(normaliseArea)
            theta = linspace(0,2.*pi,300);
            x_movingWindow = pos(1,1)+h.*cos(theta);
            y_movingWindow = pos(1,2)+h.*sin(theta);
            movingWindow = polyshape(x_movingWindow,y_movingWindow);
            AreaMovingWindow = area(intersect(ROI,movingWindow));
        end
        distances_to_other_points = Dj(i,:);
        % Removes points to far
        distances_to_other_points(i) = [];
        distances_to_other_points(distances_to_other_points>h) = [];
        lambda_j(i,1) = numel(distances_to_other_points)./AreaMovingWindow;
    end
    if(isempty(lambda_j(lambda_j(:,1)>0,1)))
        lambda_j(:) = N_j./A;
    else
        lambda_j(lambda_j(:,1)==0,1) = min(lambda_j(lambda_j(:,1)>0,1));
    end

    K = zeros(1,length(rvals_sorted));

    % Loop through distances
    for rr=1:length(rvals_sorted)
        r = rvals_sorted(rr);
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
                    m_i = m_i+winv./lambda_j(j);
                end
            end
            marks(i,rr) = (1./(lambda_i(i))).*m_i;
        end
        K(1,rr) = sum(marks(:,rr))./A;
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