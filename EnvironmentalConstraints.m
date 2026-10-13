function [newCoord,macroIdxSampled] = EnvironmentalConstraints(coord,data,macroBinLength,Rdata,Rsynthetic,...
movingType,consideredTypes,exclusionDistance,shapeData,shapeSynthetic,macroIdxSampled)

    max_try_placement = 100;

    newCoord = coord;

    warning('off', 'MATLAB:polyshape:repairedBySimplify')

    %%%%% CLUSTER RELAXATION OF MACRO-STRUCTURES
    % Step 1: Generate density map of data and of current agent placement
    [macroDensityMap_data,~,~] = buildDensityMap(data,Rdata,consideredTypes,'step',macroBinLength,'units',"counts");
    nBinMacroData = size(macroDensityMap_data,2);
    macroBinLength_exact_data = (2.*Rdata)./nBinMacroData;
    areaIntersection_data = (macroBinLength_exact_data.^2).*ones(nBinMacroData,nBinMacroData);
    if(strcmp(shapeData,"circle"))
        %%% Create mask for circle
        ROI_data_theta = linspace(0,2.*pi,360);
        ROI_data_x = Rdata.*cos(ROI_data_theta);
        ROI_data_y = Rdata.*sin(ROI_data_theta);
        ROI_data = simplify(polyshape(ROI_data_x,ROI_data_y));
        xLeft = -Rdata;
        for rr=1:nBinMacroData
            xRight = xLeft+macroBinLength_exact_data;
            yDown = -Rdata;
            for cc=1:nBinMacroData
                yUp = yDown+macroBinLength_exact_data;
                %%% Small square
                binsquare = simplify(polyshape([xLeft xLeft xRight xRight],[yUp yDown yDown yUp]));
                areaIntersection_data(rr,cc) = area(intersect(ROI_data,binsquare));
                yDown = yUp;
            end
            xLeft = xRight;
        end
    end
    insideData = areaIntersection_data>0;
    macroCountMap_data = reshape(sum(macroDensityMap_data,1),[nBinMacroData nBinMacroData]);
    macroDensityMap_data = macroCountMap_data./areaIntersection_data;
    macroDensityMap_data = macroDensityMap_data(insideData); % units: points/area unit

    [macroDensityMap_iteration,macroBinLength_exact_synthetic,agentIndexesDensityMapData] = buildDensityMap(coord,Rsynthetic,consideredTypes,'step',macroBinLength,'units',"counts");
    agentIndexesDensityMapData = agentIndexesDensityMapData(consideredTypes==movingType);
    nBinMacroSynthetic = size(macroDensityMap_iteration,2);
    areaIntersection_synthetic = (macroBinLength_exact_synthetic.^2).*ones(nBinMacroSynthetic,nBinMacroSynthetic);
    if(strcmp(shapeSynthetic,"circle"))
        %%% Create mask for circle
        ROI_synthetic_theta = linspace(0,2.*pi,360);
        ROI_synthetic_x = Rsynthetic.*cos(ROI_synthetic_theta);
        ROI_synthetic_y = Rsynthetic.*sin(ROI_synthetic_theta);
        ROI_synthetic = simplify(polyshape(ROI_synthetic_x,ROI_synthetic_y));
        xLeft = -Rsynthetic;
        for rr=1:nBinMacroSynthetic
            xRight = xLeft+macroBinLength_exact_synthetic;
            yDown = -Rsynthetic;
            for cc=1:nBinMacroSynthetic
                yUp = yDown+macroBinLength_exact_synthetic;
                %%% Small square
                binsquare = simplify(polyshape([xLeft xLeft xRight xRight],[yUp yDown yDown yUp]));
                areaIntersection_synthetic(rr,cc) = area(intersect(ROI_synthetic,binsquare));
                yDown = yUp;
            end
            xLeft = xRight;
        end
    end
    insideSynthetic = areaIntersection_synthetic>0;
    macroCountMap_iteration = reshape(sum(macroDensityMap_iteration,1),[nBinMacroSynthetic nBinMacroSynthetic]);
    macroDensityMap_iteration = macroCountMap_iteration./areaIntersection_synthetic;
    validSyntheticBins = find(insideSynthetic(:));
    validSyntheticAreas = areaIntersection_synthetic(validSyntheticBins);
    macroDensityMap_iteration_valid = macroDensityMap_iteration(validSyntheticBins);
    macroCountMap_iteration_valid = macroCountMap_iteration(validSyntheticBins);

    % Step 2: Assign desired densities to bins in current iteration
    %%% In enlarged configurations, increase density data by sampling
    increased_number = numel(macroDensityMap_iteration_valid) - numel(macroDensityMap_data);
    if(increased_number>0)
        if(isempty(macroIdxSampled))
            [~,macroIdxSampled] = datasample(macroDensityMap_data,increased_number);
        end
        macroDensityMap_data = [macroDensityMap_data; macroDensityMap_data(macroIdxSampled)];
    end
    sorted_macro_DensityMap_data = sort(macroDensityMap_data);
    % Only sort bins that re valid (i.e. inside domain)
    [sorted_macroDensityMap_iteration,idxBins] = sort(macroDensityMap_iteration_valid);
    sorted_macroCount_iteration = macroCountMap_iteration_valid(idxBins);
    sorted_validAreasSynthetic = validSyntheticAreas(idxBins);
    sorted_target_counts = sorted_macro_DensityMap_data.*sorted_validAreasSynthetic;

    % We need to make sure sum(sorted_target_counts) =
    % sum(sorted_macroCount_iteration)
    if(sum(sorted_target_counts)~=sum(sorted_macroCount_iteration))
        sorted_target_counts = sorted_target_counts./sum(sorted_target_counts);
        sorted_target_counts = ceil(sorted_target_counts.*sum(sorted_macroCount_iteration));
    end
    sortedFullBins = validSyntheticBins(idxBins);
    sortedPositionOfFullBin = zeros(nBinMacroSynthetic*nBinMacroSynthetic,1);
    sortedPositionOfFullBin(sortedFullBins) = 1:numel(sortedFullBins);

    % Step 3: Balance cells
    diffCounts = floor(sorted_macroCount_iteration-sorted_target_counts);
    overBins = find(diffCounts>0);
    visitedBins = [];
    while(~isempty(overBins))
        b = overBins(1);
        % Step 3.1: Get index of current bin
        currentBin = sortedFullBins(b);
        % Step 3.2: Check how many cells you need to redistribute
        numMovingCells = min(length(agentIndexesDensityMapData.bins(currentBin).agents),diffCounts(b));
        % Step 3.3: Choose random cells to move
        cellsInBin = agentIndexesDensityMapData.bins(currentBin).agents;
        chosen_cells = cellsInBin(randperm(length(cellsInBin), numMovingCells));
        % Step 3.4: Move cells
        for cc=chosen_cells(:).'
            % Get Moore neighborhood of current bin
            foundNeighborhood = false;
            layer = 1;
            while(~foundNeighborhood)
                neigh_bins_idx = getMooreNeighborhood(nBinMacroSynthetic,currentBin,layer);
                neigh_bins_idx = neigh_bins_idx(insideSynthetic(neigh_bins_idx));
                neigh_bins_idx_InSortedList = sortedPositionOfFullBin(neigh_bins_idx);
                diffNeigh = diffCounts(neigh_bins_idx_InSortedList);
                if(any(diffNeigh<0))
                    foundNeighborhood = true;
                else
                    layer = layer+1;
                end
            end
            % Select neighborhood bin to move to
            neighBins_underDensity_InSortedList = neigh_bins_idx_InSortedList(diffNeigh<0);
            neighBins_underDensity = neigh_bins_idx(diffNeigh<0);
            chosenBinIndex = randi(length(neighBins_underDensity));
            chosen_neigh_InSortedList = neighBins_underDensity_InSortedList(chosenBinIndex);
            chosen_neigh = neighBins_underDensity(chosenBinIndex);
            % Calculate new location for cell inside new bin
            [rI,cI] = ind2sub([nBinMacroSynthetic,nBinMacroSynthetic],chosen_neigh);
            x_left  = (rI-1)*macroBinLength_exact_synthetic - Rsynthetic;
            x_right = rI*macroBinLength_exact_synthetic - Rsynthetic;
            y_down  = (cI-1)*macroBinLength_exact_synthetic - Rsynthetic;
            y_up    = cI*macroBinLength_exact_synthetic - Rsynthetic;

            okNeighbors = false;
            placementTry = 1;
            bestProposedPosition = zeros(1,2);
            current_min_position = 10^6;
            validNeigh = ismember(newCoord(:,end), consideredTypes);
            validNeigh(cc) = false;
            neigh = newCoord(validNeigh,1:2);
            while(~okNeighbors)
                % Generate a position inside the bin
                if(strcmp(shapeSynthetic,"circle"))
                    insideChosenBin = false;
                    while(~insideChosenBin)
                        new_x = x_left + (x_right-x_left)*rand();
                        new_y = y_down + (y_up-y_down)*rand();
                        insideChosenBin = ((new_x^2 + new_y^2) <= (Rsynthetic^2));
                    end
                else
                    new_x = x_left + (x_right-x_left)*rand();
                    new_y = y_down + (y_up-y_down)*rand();
                end
                if(isempty(neigh))
                    okNeighbors = true;
                else
                    d2 = (neigh(:,1)-new_x).^2+(neigh(:,2)-new_y).^2;
                    if(min(d2)>=exclusionDistance.^2)
                        okNeighbors = true;
                    else
                        if(placementTry==1)
                            bestProposedPosition = [new_x new_y];
                            current_min_position = min(d2);
                        elseif(min(d2)>current_min_position)
                            bestProposedPosition = [new_x new_y];
                            current_min_position = min(d2);
                        end
                        placementTry = placementTry+1;
                    end
                end
                if(~okNeighbors && placementTry>=max_try_placement)
                    okNeighbors = true;
                    new_x = bestProposedPosition(1);
                    new_y = bestProposedPosition(2);
                end
            end
            newCoord(cc,1:2) = [new_x,new_y];
            % Remove cell from previous bin and add it to new bin
            sorted_macroCount_iteration(b) = sorted_macroCount_iteration(b)-1;
            sorted_macroCount_iteration(chosen_neigh_InSortedList) = sorted_macroCount_iteration(chosen_neigh_InSortedList)+1;
            sorted_macroDensityMap_iteration(b) = sorted_macroCount_iteration(b)./sorted_validAreasSynthetic(b);
            sorted_macroDensityMap_iteration(chosen_neigh_InSortedList) = sorted_macroCount_iteration(chosen_neigh_InSortedList)./sorted_validAreasSynthetic(chosen_neigh_InSortedList);
            agentIndexesDensityMapData.bins(currentBin).agents = setdiff(agentIndexesDensityMapData.bins(currentBin).agents,cc);
            agentIndexesDensityMapData.bins(chosen_neigh).agents = [agentIndexesDensityMapData.bins(chosen_neigh).agents,cc];
            % Recalculate differences
            diffCounts = floor(sorted_macroCount_iteration-sorted_target_counts);
        end
        overBins = find(diffCounts>0);
        visitedBins = [visitedBins,b];
        % Remove current bin in case condition can't be satisfied with
        % moving cell type (avoids looping forever)
        overBins = setdiff(overBins,visitedBins);
    end

end