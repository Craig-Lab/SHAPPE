function [cells,augmented_radius] = randomConfigurationGenerator(N,initialProportions,...
    areaIncrease,R,gridshape,initialShape,minDistance,voids,pAcceptCold,attractingRegions)
    if(nargin < 8)
        voids = [];
    end
    if(nargin<9)
        pAcceptCold = 0;
    end
    if(nargin<10)
        attractingRegions = [];
    end
    cellcounts = round(N.*initialProportions);
    if(sum(cellcounts)<N)
        dec = abs((cellcounts-N.*initialProportions));
        [~,ind_max_dec]=max(dec);
        cellcounts(ind_max_dec) = cellcounts(ind_max_dec)+1;
    end

    if(strcmp(initialShape,"circle"))
        initialArea = pi.*R.^2;
    elseif(strcmp(initialShape,"square"))
        initialArea = (2.*R).^2;
    end
    augmented_total_area = initialArea.*areaIncrease;
    if(strcmp(gridshape,'circle'))
        augmented_radius = sqrt(augmented_total_area./pi);
    else
        augmented_radius = sqrt(augmented_total_area)./2;
    end

    cells = zeros(N,4);
    cellOfCurrentType = 0;
    currentType = 1;
    dmin_currenttype = minDistance;
    cc=0;
    counterToFindCell = 0;
    maxTryPerCell = 1000;
    dmins = zeros(N,1);
    bestpos = zeros(1,2);
    bestdist = (-1).*10^(6);
    while(cc<N)
        if(cellOfCurrentType>=cellcounts(currentType))
            cellOfCurrentType = 0;
            currentType = currentType+1;
            dmin_currenttype = minDistance;
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
            dx = cells(1:cc,1) - x;
            dy = cells(1:cc,2) - y;
            dmin = max(dmins(1:cc),dmin_currenttype);

            okVoid = true;
            if(~isempty(voids))
                % Check if cell falls within a void
                for vv=1:size(voids,1)
                   dv = norm([x,y]-[voids(vv,1),voids(vv,2)]);
                   if(dv<voids(vv,3))
                       okVoid = false;
                        break
                   end % End of if(dv<voids(vv,3))
                end % for vv=1:size(voids,1)
            end % if(~isempty(voids))
            
            okAttractingRegions = true;
            if(~isempty(attractingRegions))
                okAttractingRegions = false;
                % Check if cell falls within an attracting regions
                for vv=1:size(attractingRegions,1)
                   dv = norm([x,y]-[attractingRegions(vv,1),attractingRegions(vv,2)]);
                   if(dv<attractingRegions(vv,3))
                       okAttractingRegions = true;
                       break
                   end % End of if(dv<attractingRegions(vv,3))
                end % for vv=1:size(attractingRegions,1)
            end % if(~isempty(attractingRegions)) ... end

            if(okAttractingRegions && okVoid)
                d2 = dx.*dx + dy.*dy;
                mind2 = min(d2);
                ok = (cc==0 | mind2>=dmin.^2);
                if(ok)
                    cells(cc+1,:) = [x,y,0,currentType];
                    cellOfCurrentType = cellOfCurrentType+1;
                    cc = cc+1;
                    counterToFindCell = 0;
                    dmins(cc,1) = dmin_currenttype;
                    bestpos = zeros(1,2);
                    bestdist = (-1).*10^(6);
                else
                    if(mind2>bestdist)
                        bestpos = [x,y];
                        bestdist = min(d2);
                    end
                    counterToFindCell = counterToFindCell+1;
                end
            else
                counterToFindCell = counterToFindCell+1;
            end
            if(counterToFindCell>=maxTryPerCell)
                x = bestpos(1);
                y = bestpos(2);
                cells(cc+1,:) = [x,y,0,currentType];
                cellOfCurrentType = cellOfCurrentType+1;
                cc = cc+1;
                counterToFindCell = 0;
                dmins(cc,1) = dmin_currenttype;
                bestpos = zeros(1,2);
                bestdist = (-1).*10^(6);
            end

        end
    end
end