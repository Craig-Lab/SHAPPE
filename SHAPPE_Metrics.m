function metrics = SHAPPE_Metrics(informationBaseline,informationComparison,rL,ripleyLweights,...
    LBAll,UBAll,minCellCount,celltypes,rneigh,metricsToCalculate)

    % To set which metrics to calculate when sin that function, enter
    % metricsToCalculate as a structure with
    % metricsToCalculate.Ldiff = true;
    % metricsToCalculate.OV_NN = true;
    % metricsToCalculate.ExposureMatrices = true;
    % metricsToCalculate.DensityMapDifference = true;

    if(nargin<10)
        calculateLdiff = true;
        calculateOV_NN = true;
        calculateExposureMatrices = true;
        calculateDensityMapDifference = true;
    else
        if(isfield(metricsToCalculate,'Ldiff'))
            calculateLdiff = metricsToCalculate.Ldiff;
        else
            calculateLdiff = true;
        end
        if(isfield(metricsToCalculate,'OV_NN'))
            calculateOV_NN = metricsToCalculate.OV_NN;
        else
            calculateOV_NN = true;
        end
        if(isfield(metricsToCalculate,'ExposureMatrices'))
            calculateExposureMatrices = metricsToCalculate.ExposureMatrices;
        else
            calculateExposureMatrices = true;
        end
        if(isfield(metricsToCalculate,'DensityMapDifference'))
            calculateDensityMapDifference = metricsToCalculate.DensityMapDifference;
        else
            calculateDensityMapDifference = true;
        end
    end

    warning('off','all')

    % Retrieve information from input structures
    %%%%% Baseline/Data
    baselineConfig = informationBaseline.configuration;
    radiusBaseline = informationBaseline.radius;
    shapeBaseline = informationBaseline.shape;
    %%%%% Comparison
    comparisonConfig = informationComparison.configuration;
    radiusComparison = informationComparison.radius;
    shapeComparison = informationComparison.shape;

    Ncelltypes = length(celltypes);

    sum_weights_celltypes = 0; % Outer loop sum of weights
    sum_weights_pair_celltypes = 0; % Inner loop sum of weights
    RipleyDifference_raw = zeros(1,length(rL));
    RipleyExcess_CI = zeros(1,length(rL));
    OV_NN = 0;
    densityMapDifference = 0;
    if(calculateLdiff || calculateOV_NN || calculateDensityMapDifference)
        if(calculateDensityMapDifference)
            %%%%%%%%%% Density Map Difference
            %%% Build density/intensity maps
            [densityMapsBaseline,~,~] = buildDensityMap(baselineConfig,radiusBaseline,celltypes,'step',80,'units',"counts");
            nbinBaseline = size(densityMapsBaseline,2);
            baselineBinLength = (2.*radiusBaseline)./nbinBaseline;
            areaIntersection_baseline = (baselineBinLength.^2).*ones(nbinBaseline,nbinBaseline);
            if(strcmp(shapeBaseline,"circle"))
                %%% Create mask for circle
                ROI_baseline_theta = linspace(0,2.*pi,360);
                ROI_baseline_x = radiusBaseline.*cos(ROI_baseline_theta);
                ROI_baseline_y = radiusBaseline.*sin(ROI_baseline_theta);
                ROI_comparison = simplify(polyshape(ROI_baseline_x,ROI_baseline_y));
                xLeft = -radiusBaseline;
                for rr=1:nbinBaseline
                    xRight = xLeft+baselineBinLength;
                    yDown = -radiusBaseline;
                    for cc=1:nbinBaseline
                        yUp = yDown+baselineBinLength;
                        %%% Small square
                        binsquare = simplify(polyshape([xLeft xLeft xRight xRight],[yUp yDown yDown yUp]));
                        areaIntersection_baseline(rr,cc) = area(intersect(ROI_comparison,binsquare));
                        yDown = yUp;
                    end
                    xLeft = xRight;
                end
            end
            insideBaseline = areaIntersection_baseline>0;
            densityMapsBaseline = reshape(sum(densityMapsBaseline,1),[nbinBaseline nbinBaseline]);
            densityMapsBaseline = densityMapsBaseline./areaIntersection_baseline;
            densityMapsBaseline = densityMapsBaseline(insideBaseline); % units: points/area unit

            [densityMapsComparison,~,~] = buildDensityMap(comparisonConfig,radiusComparison,celltypes,'step',80,'units',"counts");
            nbinComparison = size(densityMapsComparison,2);
            comparisonBinLength = (2.*radiusComparison)./nbinComparison;
            areaIntersection_comparison = (comparisonBinLength.^2).*ones(nbinComparison,nbinComparison);
            if(strcmp(shapeComparison,"circle"))
                %%% Create mask for circle
                ROI_comparison_theta = linspace(0,2.*pi,360);
                ROI_comparison_x = radiusComparison.*cos(ROI_comparison_theta);
                ROI_comparison_y = radiusComparison.*sin(ROI_comparison_theta);
                ROI_comparison = simplify(polyshape(ROI_comparison_x,ROI_comparison_y));
                xLeft = -radiusComparison;
                for rr=1:nbinComparison
                    xRight = xLeft+comparisonBinLength;
                    yDown = -radiusComparison;
                    for cc=1:nbinComparison
                        yUp = yDown+comparisonBinLength;
                        %%% Small square
                        binsquare = simplify(polyshape([xLeft xLeft xRight xRight],[yUp yDown yDown yUp]));
                        areaIntersection_comparison(rr,cc) = area(intersect(ROI_comparison,binsquare));
                        yDown = yUp;
                    end
                    xLeft = xRight;
                end
            end
            insideComparison = areaIntersection_comparison>0;
            densityMapsComparison = reshape(sum(densityMapsComparison,1),[nbinComparison nbinComparison]);
            densityMapsComparison = densityMapsComparison./areaIntersection_comparison;
            densityMapsComparison = densityMapsComparison(insideComparison); % units: points/area unit

            ovDM = overlap_ovl(densityMapsBaseline, densityMapsComparison);
            densityMapDifference = 1-ovDM;

        end % End of if(calculateDensityMapDifference)

        for nc1=1:Ncelltypes
            coordComparison = comparisonConfig(comparisonConfig(:,end)==celltypes(nc1),1:2);
            coordBaseline = baselineConfig(baselineConfig(:,end)==celltypes(nc1),1:2);
            if(size(coordBaseline,1)>=minCellCount)
                wi = ripleyLweights(nc1);
                sum_weights_celltypes = sum_weights_celltypes+wi;
                if(calculateLdiff || calculateOV_NN)
                    for nc2=1:Ncelltypes
                        Nj = length(find(baselineConfig(:,end)==celltypes(nc2)));
                        if(Nj>=minCellCount)
                            wj = ripleyLweights(nc2);
                            sum_weights_pair_celltypes = sum_weights_pair_celltypes+wi.*wj;

                            %%%%%%%%%% Ripley's H
                            if(calculateLdiff)
                                % Convert function K to H
                                LB = sqrt(max(0,reshape(LBAll(nc1,nc2,:),[1,length(rL)]))./pi)-rL;
                                UB = sqrt(reshape(UBAll(nc1,nc2,:),[1,length(rL)])./pi)-rL;
                                Hval = sqrt(reshape(informationBaseline.RipleyK(nc1,nc2,:),[1,length(rL)])./pi)-rL;
                                if(informationBaseline.InhomogeneityStatus)
                                    [L,~] = ripleyKInhom(comparisonConfig, celltypes(nc1), celltypes(nc2), radiusComparison, rL,informationBaseline.h_kernel,"H",shapeComparison);
                                else
                                    [L,~,~] = ripleyK(comparisonConfig, celltypes(nc1), celltypes(nc2), radiusComparison, rL,"H",shapeComparison);
                                end % End of if(informationBaseline.InhomogeneityStatus) ... else ... end
                                den = UB-LB;
                                den(den==0) = 1;
                                exceedance = (max(0, LB - L) + max(0, L - UB))./den;
                                RipleyExcess_CI = RipleyExcess_CI+wi.*wj.*exceedance;
                                RipleyDifference_raw = RipleyDifference_raw+wi.*wj.*abs(L-Hval);
                            end % End of if(calculateLdiff)
        
                            %%%%%%%%%% NN Overlap
                            if(calculateOV_NN)
                                NN_Baseline = getNearestNeighbors(celltypes(nc1),celltypes(nc2),baselineConfig,1);
                                NN_Comparison = getNearestNeighbors(celltypes(nc1),celltypes(nc2),comparisonConfig,1);
                                OV_NN = OV_NN+wi*wj*overlap_ovl(NN_Baseline,NN_Comparison,1024);
                            end % End of if(calculateOV_NN)
                        end % End of if(Nj>=minCellCount)
                    end % End of for nc2=1:Ncelltypes
                end % End of if(calculateLdiff || calculateOV_NN)
            end % End of if(size(coordBaseline,1)>=minCellCount)
        end % End of for nc1=1:Ncelltypes
        if(sum_weights_pair_celltypes>0 && (calculateLdiff || calculateOV_NN))
            if(calculateLdiff)
                RipleyExcess_CI = RipleyExcess_CI./sum_weights_pair_celltypes;
                RipleyDifference_raw = RipleyDifference_raw./sum_weights_pair_celltypes;
            end % End of if(calculateLdiff)
            if(calculateOV_NN)
                OV_NN = 1-OV_NN./sum_weights_pair_celltypes;
            end % End of if(calculateOV_NN)
        end % End of if(sum_weights_pair_celltypes>0 && (calculateLdiff || calculateOV_NN))
    end % End of if(calculateLdiff || calculateOV_NN || calculateDensityMapDifference)

    %%%%%%%%%% Exposure matrices
    score_mu = 0;
    score_pi = 0;
    score_e = 0;
    if(calculateExposureMatrices)        
        S_Baseline   = neighborhoodExposureStats(baselineConfig(:,[1 2 end]),  celltypes, rneigh);
        S_Comparison  = neighborhoodExposureStats(comparisonConfig(:,[1 2 end]), celltypes, rneigh);
        for t=1:length(celltypes)
            Ni = length(find(baselineConfig(:,end)==celltypes(t)));
            if(Ni<minCellCount)
                S_Baseline.mu(t,:) = nan;
                S_Baseline.mu(:,t) = nan;
                S_Baseline.pi(t,:) = nan;
                S_Baseline.pi(:,t) = nan;
                S_Baseline.e(t,:) = nan;
                S_Baseline.e(:,t) = nan;
                S_Comparison.mu(t,:) = nan;
                S_Comparison.mu(:,t) = nan;
                S_Comparison.pi(t,:) = nan;
                S_Comparison.pi(:,t) = nan;
                S_Comparison.e(t,:) = nan;
                S_Comparison.e(:,t) = nan;
            end % if(Ni<minCellCount)
        end % End of for t=1:length(celltypes)
        score_mu  = sum(abs(S_Baseline.mu(:) - S_Comparison.mu(:)),"omitnan");
        score_pi  = sum(abs(S_Baseline.pi(:) - S_Comparison.pi(:)),"omitnan");
        score_e   = sum(abs(S_Baseline.e(:)  - S_Comparison.e(:)),"omitnan");
    end % End of if(calculateExposureMatrices)


    metrics = struct();
    metrics.RipleyExcess_CI = RipleyExcess_CI;
    metrics.RipleyDifference_raw = RipleyDifference_raw;
    metrics.OV_NN = OV_NN;
    metrics.densityMapDifference = densityMapDifference;
    metrics.score_mu = score_mu;
    metrics.score_pi = score_pi;
    metrics.score_e = score_e;

end