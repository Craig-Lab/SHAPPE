function distancesConfiguration = getDistanceDistributions(configurationData,...
    pointtypeIndexes,pointtypeNames,numNeigh)

    if(nargin<4)
        numNeigh = 1;
    end

    Ntypes = numel(pointtypeIndexes);
    distancesConfiguration = struct();

    % Outer loop
    for t1=1:Ntypes
        ind1 = pointtypeIndexes(t1);
        type1 = pointtypeNames(t1);
        % Inner loop
        for t2=1:Ntypes
            ind2 = pointtypeIndexes(t2);
            type2 = pointtypeNames(t2);
            current_dist = getNearestNeighbors(ind1,ind2,configurationData,numNeigh);
            if(~isempty(current_dist))
                if(numNeigh == 1 && iscolumn(current_dist))
                    current_dist = current_dist.'; 
                end
            end
            if(numNeigh>1 && ~isempty(current_dist))
                current_dist = [current_dist(1,:);diff(current_dist)];
            end
            distancesConfiguration.(type1).(type2).data = current_dist;
            for nr=1:numNeigh
                d = current_dist(nr,:);
                d = d(~isnan(d) & ~isinf(d));
            
                if isempty(d)
                    distancesConfiguration.(type1).(type2).distribution{nr} = [];
                    distancesConfiguration.(type1).(type2).cum_probs{nr}    = [];
                    continue
                end
            
                % Sort
                d = sort(d);
            
                % Build CDF
                [xU, ~, ic] = unique(d, 'stable');
                counts = accumarray(ic, 1);
                F = cumsum(counts)./numel(d);
            
                % Store
                distancesConfiguration.(type1).(type2).distribution{nr} = xU(:).';
                distancesConfiguration.(type1).(type2).cum_probs{nr}    = F(:).';
            end % End of for nr=1:numNeigh

        end % End of for t2=1:length(pointtypeIndexes)
        distancesConfiguration.(type1).center = (sqrt((configurationData(configurationData(:,end)==ind1,1).^2)+(configurationData(configurationData(:,end)==ind1,2).^2)))';
        distancesConfiguration.(type1).x = (configurationData(configurationData(:,end)==ind1,1))';
        distancesConfiguration.(type1).y = (configurationData(configurationData(:,end)==ind1,2))';
        distancesConfiguration.(type1).R = sqrt(((configurationData(configurationData(:,end)==ind1,1))').^2+((configurationData(configurationData(:,end)==ind1,2))').^2);
        distancesConfiguration.(type1).theta = atan2((configurationData(configurationData(:,end)==ind1,2))',(configurationData(configurationData(:,end)==ind1,1))');
    end % End of for t1=1:length(pointtypeIndexes)
   
end

