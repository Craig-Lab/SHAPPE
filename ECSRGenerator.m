function [points,augmented_radius,repulsiveRegions,attractiveRegions,pointcounts] = ...
    ECSRGenerator(N,initialProportions,areaIncrease,RObs,gridshape,observedShape,...
    exclusionDistance,types,homogeneityStatus,ECSRAttractionScenario,LowDensityRegionsAnalysis)

    %%%%% STEP 0: Check you have enough information
    if(~homogeneityStatus && nargin<11)
        msg = "You must provide information for ECSRAttractionScenario and LowDensityRegionAnalysis for ECSR.";
        error(msg)
    end

    %%%%% STEP 1: Calculate new point counts and radius
    rawCounts = N.*initialProportions;
    pointcounts = floor(rawCounts);
    [~,order] = sort(rawCounts - pointcounts, 'descend');
    k = N - sum(pointcounts);
    pointcounts(order(1:k)) = pointcounts(order(1:k)) + 1;

    if(strcmp(observedShape,"circle"))
        initialArea = pi.*RObs.^2;
    elseif(strcmp(observedShape,"square"))
        initialArea = (2.*RObs).^2;
    end
    augmented_total_area = initialArea.*areaIncrease;
    if(strcmp(gridshape,'circle'))
        augmented_radius = sqrt(augmented_total_area./pi);
    else
        augmented_radius = sqrt(augmented_total_area)./2;
    end

    %%%%% STEP 2: Generate repulsive/attractive regions if inhomogeneous
    if(homogeneityStatus)
        repulsiveRegions = [];
        attractiveRegions =[];
    else
        max_tries = 100;
        minDistanceRegions = 20;
        if(~ECSRAttractionScenario)
            attractiveRegions = [];
            coldAreas = LowDensityRegionsAnalysis.AreaDistribution;
            NcoldObs = numel(coldAreas);
            Ncold = round(NcoldObs.*areaIncrease);
            % Generate repulsive regions
            repulsiveRegions = zeros(Ncold,3);
            v = 1;
            current_try = 1;
            while(v<=Ncold)
                % Sample effective radius
                if(v<=NcoldObs)
                    area = coldAreas(v);
                else
                    area = coldAreas(randi(NcoldObs));
                end
                r = sqrt(area./pi);
                % Sample center
                thetaPos = 2*pi*rand();
                rPos = (augmented_radius-r)*sqrt(rand());
                x = rPos*cos(thetaPos);
                y = rPos*sin(thetaPos);
                if(v==1)
                    repulsiveRegions(v,:) = [x,y,r];
                    v = v+1;
                else
                    % Check if overlapping with already placed voids
                    ok = true;
                    for vv=1:(v-1)
                        x_other = repulsiveRegions(vv,1);
                        y_other = repulsiveRegions(vv,2);
                        r_other = repulsiveRegions(vv,3);
                        d = norm([x,y]-[x_other,y_other]);
                        if(d<=(r_other+r+minDistanceRegions))
                            ok = false;
                            break
                        end % End of if(d<=(r_other+r))
                    end % End of for vv=1:(v-1)
                    if(ok)
                        repulsiveRegions(v,:) = [x,y,r];
                        v = v+1;
                        current_try = 1;
                    else
                        current_try = current_try+1;
                        if(current_try>max_tries)
                            current_try = 1;
                            v = 1;
                        end % End of if(current_try>max_tries)
                    end % End of if(ok)
                end % End of if(v==1)
            end % End of while(v<=Ncold)
        else
            repulsiveRegions = [];
            hotAreas = LowDensityRegionsAnalysis.AreaDistributionHot;
            NhotObs = numel(hotAreas);
            Nhot = round(NhotObs.*areaIncrease);
            % Generate attracting regions
            attractiveRegions = zeros(Nhot,3);
            v = 1;
            current_try = 1;
            while(v<=Nhot)
                % Sample effective radius
                if(v<=NhotObs)
                    area = hotAreas(v);
                else
                    area = hotAreas(randi(NhotObs));
                end
                r = sqrt(area./pi);
                % Sample center
                if(r<augmented_radius)
                    thetaPos = 2*pi*rand();
                    rPos = (augmented_radius-r)*sqrt(rand());
                    x = rPos*cos(thetaPos);
                    y = rPos*sin(thetaPos);
                else
                    x = 0;
                    y = 0;
                end
                if(v==1)
                    attractiveRegions(v,:) = [x,y,r];
                    v = v+1;
                else
                    % Check if overlapping with already placed
                    % attracting regions
                    ok = true;
                    for vv=1:(v-1)
                        x_other = attractiveRegions(vv,1);
                        y_other = attractiveRegions(vv,2);
                        r_other = attractiveRegions(vv,3);
                        d = norm([x,y]-[x_other,y_other]);
                        if(d<=(r_other+r+minDistanceRegions))
                            ok = false;
                            break
                        end % End of if(d<=(r_other+r))
                    end % End of for vv=1:(v-1)
                    if(ok)
                        attractiveRegions(v,:) = [x,y,r];
                        v = v+1;
                        current_try = 1;
                    else
                        current_try = current_try+1;
                        if(current_try>max_tries)
                            current_try = 1;
                            v = 1;
                        end % End of if(current_try>max_tries)
                    end % End of if(ok)
                end % End of if(v==1)
            end % End of while(v<Nhot)
        end % end of if(~ECSRAttractionScenario) ... else ... end
    end % end of if(homogeneityStatus) ... else ... end

    %%%%% STEP 3: Place points
    points = zeros(N,4);
    pointsOfCurrentType = 0;
    currentType = 1;
    dmin_currenttype = exclusionDistance;
    cc=0;
    counterToFindPosition = 0;
    maxTryPerPoint = 1000;
    dmins = zeros(N,1);
    bestpos = zeros(1,2);
    bestdist = (-1).*10^(6);
    while(cc<N)
        if(pointsOfCurrentType>=pointcounts(currentType))
            pointsOfCurrentType = 0;
            currentType = currentType+1;
            dmin_currenttype = exclusionDistance;
        else
            if(strcmp(gridshape,'circle'))
                theta = 2*pi*rand();
                r = augmented_radius*sqrt(rand());
                x = r*cos(theta);
                y = r*sin(theta);
            else
                x = unifrnd(-augmented_radius,augmented_radius);
                y = unifrnd(-augmented_radius,augmented_radius);
            end
            dx = points(1:cc,1) - x;
            dy = points(1:cc,2) - y;
            dmin = max(dmins(1:cc),dmin_currenttype);

            okRepulsiveRegions = true;
            if(~isempty(repulsiveRegions))
                % Check if cell falls within a void
                for vv=1:size(repulsiveRegions,1)
                   dv = norm([x,y]-[repulsiveRegions(vv,1),repulsiveRegions(vv,2)]);
                   if(dv<repulsiveRegions(vv,3))
                       okRepulsiveRegions = false;
                        break
                   end % End of if(dv<repulsiveRegions(vv,3))
                end % for vv=1:size(repulsiveRegions,1)
            end % if(~isempty(repulsiveRegions))
            
            okAttractiveRegions = true;
            if(~isempty(attractiveRegions))
                okAttractiveRegions = false;
                % Check if cell falls within an attractive regions
                for vv=1:size(attractiveRegions,1)
                   dv = norm([x,y]-[attractiveRegions(vv,1),attractiveRegions(vv,2)]);
                   if(dv<attractiveRegions(vv,3))
                       okAttractiveRegions = true;
                       break
                   end % End of if(dv<attractiveRegions(vv,3))
                end % for vv=1:size(attractiveRegions,1)
            end % if(~isempty(attractiveRegions)) ... end

            if(okAttractiveRegions && okRepulsiveRegions)
                d2 = dx.*dx + dy.*dy;
                mind2 = min(d2);
                ok = (cc==0 || all(d2>=dmin.^2));
                if(ok)
                    points(cc+1,:) = [x,y,0,types(currentType)];
                    pointsOfCurrentType = pointsOfCurrentType+1;
                    cc = cc+1;
                    counterToFindPosition = 0;
                    dmins(cc,1) = dmin_currenttype;
                    bestpos = zeros(1,2);
                    bestdist = (-1).*10^(6);
                else
                    if(mind2>bestdist)
                        bestpos = [x,y];
                        bestdist = min(d2);
                    end
                    counterToFindPosition = counterToFindPosition+1;
                end
            else
                counterToFindPosition = counterToFindPosition+1;
            end
            if(counterToFindPosition>=maxTryPerPoint)
                x = bestpos(1);
                y = bestpos(2);
                points(cc+1,:) = [x,y,0,types(currentType)];
                pointsOfCurrentType = pointsOfCurrentType+1;
                cc = cc+1;
                counterToFindPosition = 0;
                dmins(cc,1) = dmin_currenttype;
                bestpos = zeros(1,2);
                bestdist = (-1).*10^(6);
            end

        end
    end
end