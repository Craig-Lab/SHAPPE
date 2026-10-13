function densityGridAnalysis = densityGridAnalysis(conf,R,binLength,shape,thresholdInformation)

    DensityThreshold = thresholdInformation.DensityThreshold;

    if(isfield(thresholdInformation,'minArea'))
        minArea = thresholdInformation.minArea;
    else
        minArea = pi*(50).^2;
    end

    nBins = round(2.*R./binLength);
    bL = 2.*R./nBins;

    Ncells = size(conf,1);
    areaConf = (2.*R).^2;
    dens0 = Ncells./areaConf;
    if(isfield(thresholdInformation,'quantile'))
        DensityThreshold = dens0.*thresholdInformation.quantile;
    end % End of if(isfield(thresholdInformation,'quantile'))

    if(strcmp(shape,"circle"))
        % Calculate intersection between a unit circle and the grid inscribed onto
        % it with nBin = nbinIMC
        theta = linspace(0, 2*pi, 360);
        x_unitcircle = cos(theta);
        y_unitcircle = sin(theta);
        unitCircle = polyshape(x_unitcircle, y_unitcircle);

        propIntersectionsConf = zeros(nBins,nBins);
        bLunit = 2/nBins;
        y_up = 1;
        for rr=1:nBins
            x_left = -1;
            y_down = y_up-bLunit;
            for cc=1:nBins
                x_right = x_left+bLunit;
                current_square = polyshape([x_left x_left x_right x_right],[y_up y_down y_down y_up]);
                intersectionSquareCircle = intersect(unitCircle, current_square);
                area_intersection = area(intersectionSquareCircle);
                propIntersectionsConf(rr,cc) = area_intersection./(bLunit.^2);
                x_left = x_right;
            end % end of for cc=1:nBinsConf1
            y_up = y_down;
        end % end of for rr=1:nBinsConf1
    elseif(strcmp(shape,"square"))
        propIntersectionsConf = ones(nBins,nBins);
    end

    x = conf(:,1);
    y = conf(:,2);

    D = zeros(nBins,nBins);
    y_up = R;
    for rr=1:nBins
        x_left = -R;
        y_down = y_up-bL;
        for cc=1:nBins
            x_right = x_left+bL;
            inX = (x >= x_left) & ( (cc < nBins & x < x_right) | (cc == nBins & x <= x_right) );
            inY = (y >= y_down) & ( (rr < nBins & y < y_up)    | (rr == nBins & y <= y_up)    );
            num_cells_inside = sum(inX & inY);
            D(rr,cc) = num_cells_inside./(bL.*bL.*propIntersectionsConf(rr,cc));
            x_left = x_right;
        end % end of for cc=1:nBinsConf1
        y_up = y_down;
    end % end of for rr=1:nBinsConf1

    % Analyse the number and areas of cold sections
    min_number_connected = ceil(minArea./(bL.^2));
    binary_mask_UnderDensityThreshold = (D<=DensityThreshold);
    AreaDistributionCold = [];
    conn = 8;
    BW_raw = binary_mask_UnderDensityThreshold;
    BW = binary_mask_UnderDensityThreshold;
    if(min_number_connected>1)
        % XXX Change this fo circle
        BW = bwareaopen(BW, min_number_connected, conn);
    end
    CC = bwconncomp(BW,conn);
    num_patches_all_included = CC.NumObjects;

    hot_regions_maks = bwareaopen(~BW & propIntersectionsConf>0, min_number_connected, conn);
    CC_hot = bwconncomp(hot_regions_maks,conn);
    num_hot_patches = CC_hot.NumObjects;
    AreaDistributionHot = [];

    largest_UnderDensityThreshold_region_area = 0;
    average_UnderDensityThreshold_region_area = 0;
    num_patches = 0;
    if num_patches_all_included > 0
        patchAreas = nan(num_patches_all_included,1);
        for i = 1:num_patches_all_included
            idx = CC.PixelIdxList{i}; % linear indices into the mask/grid
            current_area = sum(propIntersectionsConf(idx)) * (bL*bL);
            % Sum true area of bins in this patch (handles circle intersections too)
            if(current_area>=minArea)
                patchAreas(i) = sum(propIntersectionsConf(idx)) * (bL*bL);
                AreaDistributionCold = [AreaDistributionCold,patchAreas(i)];
                num_patches = num_patches+1;
            end
        end
        valid = ~isnan(patchAreas);
        if any(valid)
            largest_UnderDensityThreshold_region_area = max(patchAreas(valid));
            average_UnderDensityThreshold_region_area = mean(patchAreas(valid));
        end
    end

    if num_hot_patches > 0
        for i = 1:num_hot_patches
            idx = CC_hot.PixelIdxList{i}; % linear indices into the mask/grid
            current_area = sum(propIntersectionsConf(idx)) * (bL*bL);
            % Sum true area of bins in this patch (handles circle intersections too)
            if(current_area>=minArea)
                AreaDistributionHot = [AreaDistributionHot,current_area];
            end
        end
    end

    densityGridAnalysis = struct();
    densityGridAnalysis.unsmoothDensityMap = D;
    densityGridAnalysis.binary_mask_UnderDensityThreshold = BW;
    densityGridAnalysis.binary_mask_UnderDensityThreshold_raw = BW_raw;
    densityGridAnalysis.CC = CC;
    densityGridAnalysis.num_patches = num_patches;
    densityGridAnalysis.average_UnderDensityThreshold_region_area = average_UnderDensityThreshold_region_area;
    densityGridAnalysis.largest_UnderDensityThreshold_region_area = largest_UnderDensityThreshold_region_area;
    densityGridAnalysis.threshold = DensityThreshold;
    densityGridAnalysis.nBins = nBins;
    densityGridAnalysis.binLength = bL;
    densityGridAnalysis.propIntersectionsConf = propIntersectionsConf;
    densityGridAnalysis.AreaDistribution = AreaDistributionCold;
    densityGridAnalysis.AreaDistributionHot = AreaDistributionHot;

end