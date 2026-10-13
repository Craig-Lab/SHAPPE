function [pointPlacement,cell_diplacement_med,endCriterionAchieved,net_displacements_med,counter,macroIdxSampled] = ...
    IterativePointPlacement(initialConfig,Xobs,targetDistances,allTypes,movingType,compare_with_types,comparisonWeights,...
    placedTypes,Rdata,Rsynthetic,shapeData,shapeSynthetic,varargin)

    % Parse through inputs
    inputStruct = inputParser;
    addRequired(inputStruct,'initialConfig');
    addRequired(inputStruct,'Xobs');
    addRequired(inputStruct,'targetDistances');
    addRequired(inputStruct,'allTypes');
    addRequired(inputStruct,'movingType');
    addRequired(inputStruct,'compare_with_types');
    addRequired(inputStruct,'comparisonWeights');
    addRequired(inputStruct,'placedTypes');
    addRequired(inputStruct,'Rdata');
    addRequired(inputStruct,'Rsynthetic');
    addRequired(inputStruct,'shapeData');
    addRequired(inputStruct,'shapeSynthetic');

    addParameter(inputStruct,'maxIterations',10000);
    addParameter(inputStruct,'phi',10);
    addParameter(inputStruct,'phiMinimum',10);
    addParameter(inputStruct,'rW',10);
    addParameter(inputStruct,'L',1);
    addParameter(inputStruct,'a',0);
    addParameter(inputStruct,'repulsiveRegions',[]);
    addParameter(inputStruct,'attractiveRegions',[]);
    addParameter(inputStruct,'macroIdxSampled',[]);
    addParameter(inputStruct,'macroBinLength',50);
    addParameter(inputStruct,'lambda',1);
    addParameter(inputStruct,'ghostDistance',inf);
    addParameter(inputStruct,'exclusionDistance',0);
    addParameter(inputStruct,'Nenv',100);
    addParameter(inputStruct,'dMed',0.0175);
    addParameter(inputStruct,'d90',10.75);
    addParameter(inputStruct,'dNet',10.75);
    addParameter(inputStruct,'nEq',210);

    parse(inputStruct, initialConfig,Xobs,targetDistances,allTypes,movingType,compare_with_types,comparisonWeights,...
    placedTypes,Rdata,Rsynthetic,shapeData,shapeSynthetic,varargin{:});

    maxIterations = inputStruct.Results.maxIterations;
    phi = inputStruct.Results.phi;
    phiMinimum = inputStruct.Results.phiMinimum;
    rW = inputStruct.Results.rW;
    L = inputStruct.Results.L;
    a = inputStruct.Results.a;
    repulsiveRegions = inputStruct.Results.repulsiveRegions;
    attractiveRegions = inputStruct.Results.attractiveRegions;
    macroIdxSampled = inputStruct.Results.macroIdxSampled;
    macroBinLength = inputStruct.Results.macroBinLength;
    lambda = inputStruct.Results.lambda;
    ghostDistance = inputStruct.Results.ghostDistance;
    exclusionDistance = inputStruct.Results.exclusionDistance;
    Nenv = inputStruct.Results.Nenv;
    dMed = inputStruct.Results.dMed;
    d90 = inputStruct.Results.d90;
    dNet = inputStruct.Results.dNet;
    nEq = inputStruct.Results.nEq;


    if(~(strcmp(shapeSynthetic,"square") || strcmp(shapeSynthetic,"circle")) || ~(strcmp(shapeData,"square") || strcmp(shapeData,"circle")))
        errMsg = "The algorithm as it is now designed only supports 'square' or 'circle' for inputs shapeData and shapeSynthetic.";
        error(errMsg)
    end 

    pointPlacement = initialConfig;

    N_agents = size(pointPlacement,1);

    % Find agents to move around
    index_agents_to_place = find(pointPlacement(:,end)==movingType);

    continueTime = true;
    counter = 0;

    % Keep track of net displacement over last number_iterations_net
    % iterations
    number_iterations_net = 25;
    displacements_net_iterations = zeros(N_agents,2,number_iterations_net);

    while(continueTime)

        displacements = zeros(N_agents,2);

        % Jacobi update
        pointPlacement0 = pointPlacement;
        %%% Positions
        XY  = pointPlacement0(:,1:2);
        %%% Agent types
        T   = pointPlacement0(:,end);

        % Create ghosts points
        symPoints = getSymmetricPoints(pointPlacement0, Rsynthetic, ghostDistance,shapeSynthetic);
        if(strcmp(shapeSynthetic,"square"))
            symPoints = symPoints( abs(symPoints(:,1)) > Rsynthetic | abs(symPoints(:,2)) > Rsynthetic , :);
        elseif(strcmp(shapeSynthetic,"circle"))
            symPoints = symPoints( symPoints(:,1).^2+symPoints(:,2).^2 > Rsynthetic.^2, :);
        end
        % XYall and Tall include agents and ghost agents
        XYall = [XY; symPoints(:,1:2)];
        Tall = [T;  symPoints(:,end)];
        
        % Create an array to store proposed new positions
        newXY = XY;

        % Create tree to search for neighbors
        Ntypes = numel(allTypes);
        trees   = cell(1, Ntypes);
        XYtype  = cell(1, Ntypes);
        idxType = cell(1, Ntypes); 
        for tt = 1:Ntypes
            tcode = allTypes(tt);
            idxAll = find(Tall == tcode);
            idxType{tt} = idxAll;
            if(~isempty(idxAll))
                XYtype{tt} = XYall(idxAll,:);
                trees{tt}  = KDTreeSearcher(XYtype{tt});
            else
                XYtype{tt} = [];
                trees{tt}  = [];
            end % End of if ~isempty(idxAll) ... else ... end
        end % End of for tt = 1:Ntypes

        % Compute update for every agent
        for nC=(index_agents_to_place)'

            displacement_vector = zeros(1,2);
            % Position of current agent
            x = pointPlacement0(nC,1);
            y = pointPlacement0(nC,2);
            % Position of nearest numNeighbors agents of every type
            nearest_agents = nan(Ntypes, L, 2);

            for tt = 1:Ntypes
                tcode = allTypes(tt);
                if(isempty(XYtype{tt}))
                    continue;
                end % End of if(isempty(XYtype{tt}))
            
                kk = min(L + 1, size(XYtype{tt},1));
                % Drop self
                nnLocal = knnsearch(trees{tt}, XY(nC,:), 'K', kk);
                isSelf  = idxType{tt}(nnLocal) == nC;
                pos = XYtype{tt}(nnLocal(~isSelf),:);
                % In case numNeighbors > size(pos,1)
                kk2 = min(L, size(pos,1));
                if(kk2 > 0)
                    nearest_agents(tt, 1:kk2, :) = pos(1:kk2,:);
                end % End of if(kk2 > 0)
            end % End of for tt = 1:Ntypes

            agent_displacement = zeros(1,2);

            % Micro-structure: Cell-cell interactions
            for nT=compare_with_types
                nT_idx = find(allTypes==nT);
                for k=1:L
                    pos_neigh = [nearest_agents(nT_idx,k,1) nearest_agents(nT_idx,k,2)];
                    if(all(~isnan(pos_neigh)))
                        displacement_vector(1,1) = pos_neigh(1)-x;
                        displacement_vector(1,2) = pos_neigh(2)-y;

                        current_min = norm(displacement_vector);
                        desired_min_distance = targetDistances(nC,L*(nT_idx-1)+k);

                        if(current_min>=exclusionDistance)
                            agent_displacement = agent_displacement+...
                                (1./(k.^a)).*comparisonWeights(nT_idx).*displacement_vector.*...
                                (current_min-desired_min_distance)./(current_min+eps);
                        elseif(current_min < 1e-12)
                            u = randn(1,2);
                            u = u/(norm(u)+eps);
                            agent_displacement = agent_displacement - rW.*u.*(exclusionDistance - current_min);
                        else
                            displacement_vector = -displacement_vector;
                            agent_displacement = agent_displacement + rW.*displacement_vector.*...
                                (exclusionDistance - current_min)./(current_min+eps);
                        end % End of if(current_min>=deltaMin) ... elseif(current_min < 1e-12) ... else ... end
                    end % End of if(all(~isnan(pos_neigh)))
                end % End of for k=1:numNeighbors
            end % End of for nT=current_compare_with

            agent_displacement = lambda.*agent_displacement;

            % Macro-structure: Environmental factors defined through
            % attractive regions
            if(~isempty(attractiveRegions))
                okAttractiveRegions = false;
                closest_attracting_region = 0;
                min_dist_to_boundary = inf;
                for vv=1:size(attractiveRegions,1)
                    % Distance to attractive region center
                    dv = norm([x,y]-[attractiveRegions(vv,1),attractiveRegions(vv,2)]);            
                    % Check if already inside an attractive region
                    if(dv<attractiveRegions(vv,3))
                        okAttractiveRegions = true;
                        break
                    else
                        % Distance to boundary of attractive region
                        dist_to_boundary = dv-attractiveRegions(vv,3);
                        if(dist_to_boundary<min_dist_to_boundary)
                            min_dist_to_boundary = dist_to_boundary; % How far from boundary
                            closest_attracting_region = vv;
                            dv_closest = dv; % How far from center
                        end % End of if(dist_to_boundary<min_dist_to_boundary)
                    end % End of if(dv<attractive(vv,3)) ... else ... end
                end % End of for vv=1:size(attractive,1)
            
                if(~okAttractiveRegions)
                    attraction_dir = [attractiveRegions(closest_attracting_region,1),attractiveRegions(closest_attracting_region,2)] - [x,y];
                    attraction_dir = attraction_dir./norm(attraction_dir);
                    agent_displacement = agent_displacement+rW.*attraction_dir.*(1 - min_dist_to_boundary./dv_closest);
                end % End of if(~okAttractingRegions)
            end % End of if(~isempty(attractionRegions))

            % Macro-structure: Environmental factors defined through
            % repulsive regions
            if(~isempty(repulsiveRegions))
                for vv=1:size(repulsiveRegions,1)
                    % Check if inside repulsiveRegions
                    dv = norm([x,y]-[repulsiveRegions(vv,1),repulsiveRegions(vv,2)]);
                    if(dv<repulsiveRegions(vv,3))
                        % Add a repulsion
                        if(dv>0)
                            void_repulsion_dir = [x,y]-[repulsiveRegions(vv,1),repulsiveRegions(vv,2)];
                        else
                            theta = unifrnd(0,2.*pi);
                            void_repulsion_dir = [repulsiveRegions(vv,3).*cos(theta) repulsiveRegions(vv,3).*sin(theta)];
                        end
                        void_repulsion_dir = void_repulsion_dir./norm(void_repulsion_dir);
                        agent_displacement = agent_displacement+rW.*void_repulsion_dir.*(1-(dv./repulsiveRegions(vv,3)));
                        break
                    end % End of if(dv<repulsiveRegions(vv,3))
                end % End of for vv=1:size(repulsiveRegions,1)
            end % End of if(~isempty(repulsiveRegions))

            if(phi>0)
                step_limit = max(phiMinimum, phi.*(1 - counter/maxIterations));
                current_size_displacement = norm(agent_displacement);
                if current_size_displacement > step_limit
                    agent_displacement = (step_limit./current_size_displacement).*agent_displacement;
                end
            end

            temp_x = agent_displacement(1,1)+x;
            temp_y = agent_displacement(1,2)+y;

            if(strcmp(shapeSynthetic,"square"))
                if(temp_x>=Rsynthetic)
                   temp_x = temp_x-2.*Rsynthetic;
                elseif(temp_x<=-Rsynthetic)
                   temp_x = temp_x+2.*Rsynthetic;
                end % End of if(temp_x>=Rsynthetic) ... elseif(temp_x<=-Rsynthetic) ... end
                if(temp_y>=Rsynthetic)
                   temp_y = temp_y-2.*Rsynthetic;
                elseif(temp_y<=-Rsynthetic)
                   temp_y = temp_y+2.*Rsynthetic;
                end % End of if(temp_y>=Rsynthetic) ... elseif(temp_y<=-Rsynthetic) ... end
                agent_displacement(1,1) = temp_x-x;
                agent_displacement(1,2) = temp_y-y;
            else
                p = [temp_x temp_y];
                r = norm(p);
                if(r>Rsynthetic)
                    n = p./r;
                    temp_x = temp_x-2.*Rsynthetic.*n(1);
                    temp_y = temp_y-2.*Rsynthetic.*n(2);
                    agent_displacement(1,1) = temp_x-x;
                    agent_displacement(1,2) = temp_y-y;
                end % End of if(Rsynthetic>Rsynthetic)
            end

            new_pos = [temp_x temp_y];
            newXY(nC,:) = new_pos;
            displacements(nC,:) = agent_displacement;
        end

        % Apply all updates
        pointPlacement(:,1:2) = newXY;

        % Update net displacements
        if (counter+1) <= number_iterations_net
            displacements_net_iterations(:,:,counter+1) = displacements;
        else
            displacements_net_iterations(:,:,1:end-1) = displacements_net_iterations(:,:,2:end);
            displacements_net_iterations(:,:,end)     = displacements;
        end

        if( counter>0 && mod(counter,Nenv)==0 && ~isempty(macroBinLength))
            prevPointPlacement = pointPlacement;
            [pointPlacement,macroIdxSampled] = EnvironmentalConstraints(prevPointPlacement,Xobs,macroBinLength,Rdata,Rsynthetic,...
                movingType,[placedTypes,movingType],exclusionDistance,shapeData,shapeSynthetic,macroIdxSampled);
            displacements_net_iterations(:,:,end) = pointPlacement(:,1:2)-prevPointPlacement(:,1:2);
        else
            displacements_size = sqrt(displacements(index_agents_to_place,1).^2 + displacements(index_agents_to_place,2).^2);
            cell_diplacement_med = median(displacements_size);
            cell_diplacement_Q90 = quantile(displacements_size,0.9);

            net_displacements = vecnorm(sum(displacements_net_iterations(index_agents_to_place,:,:),3),2,2);
            net_displacements_med = median(net_displacements);

            if(cell_diplacement_med<=dMed && cell_diplacement_Q90 <= d90)
                endCriterionAchieved = 1;
                continueTime = false;
            elseif(net_displacements_med<=dNet && (counter+1)>=nEq)
                endCriterionAchieved = 2;
                continueTime = false;
            end
            if(counter>=maxIterations)
                endCriterionAchieved = 3;
                continueTime = false;
            end
        end

        counter = counter+1;

        if(any(any(isnan(pointPlacement))))
            idx_nan_points = any(isnan(pointPlacement),2);
            disp("Points with NANs are "+...
                [sprintf('%d,', idx_nan_points(1:end-1)), sprintf('%d', idx_nan_points(end))]+"."+newline)
            data_to_print = targetDistances(idx_nan_points,:);
            for i = 1:size(data_to_print,1)
                formatSpec = repmat('%d ', 1, size(data_to_print,2));
                formatSpec = [formatSpec(1:end-1), '\n\n'];
                fprintf(formatSpec, data_to_print(i, :));
            end
            msgError = "There was a NaN in the point placement.";
            error(msgError)
        end
    end
end