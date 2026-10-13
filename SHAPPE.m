function [outputs] = SHAPPE(Xobs,options)

    % First thing first: inputs
    arguments
        Xobs (:,:) double {mustBeNonempty}

        options.RObs    (1,1) double = -1
        options.AreaObs (1,1) double = -1
        options.ShapeObs (1,1) string = "square"

        options.exclusionDistance (1,1) double {mustBeNonnegative} = 0
        options.Nmin (1,1) double {mustBeInteger,mustBePositive} = 15
        options.types double = []
        options.types_name string = string.empty
        options.typeOrdering = "descend"

        options.L (1,1) double {mustBeInteger,mustBePositive} = 1
        options.a (1,1) double {mustBeNonnegative} = 0.5
        options.rW (1,1) double {mustBeNonnegative} = 10
        options.lambda (1,1) double {mustBePositive} = 0.01/0.133
        options.phi (1,1) double {mustBeNonnegative} = 2*10.75
        options.phiMinimum (1,1) double {mustBeNonnegative} = 10.75
        options.ghostDistance (1,1) double {mustBePositive} = 50
        options.C (1,1) string = "gaussian"
        options.wNoCorr (1,1) double {mustBeNonnegative} = 0

        options.dMed (1,1) double {mustBePositive} = 0.1075
        options.d90  (1,1) double {mustBePositive} = 10.75
        options.dNet (1,1) double {mustBePositive} = 10.75
        options.nEq  (1,1) double {mustBeInteger,mustBePositive} = 210
        options.MaxN (1,1) double {mustBeInteger,mustBePositive} = 10000

        options.nEnv (1,1) double {mustBeInteger,mustBePositive} = 50
        options.bLEnv (1,1) double {mustBePositive} = 80
        options.useRegions (1,1) logical = true
        options.attract (1,1) logical = false
        options.SHAPPEalpha (1,1) double {mustBePositive,mustBeLessThan(options.SHAPPEalpha,1)} = 0.05

        options.N (1,1) double {mustBeInteger,mustBePositive} = 10000
        options.Nreps (1,1) double {mustBeInteger,mustBePositive} = 1
        options.shapeSynthetic (1,1) string = "circle"
        options.parallelLoops (1,1) double {mustBeInteger,mustBeNonnegative} = 0
        options.outDirectoryRandomBaseline string = string.empty
        options.outDirectorySHAPPE string = string.empty
        options.outDirectoryRepulsiveRegions string = string.empty
        options.outDirectoryAttractiveRegions string = string.empty

        options.alphaFtest (1,1) double {mustBePositive,mustBeLessThan(options.alphaFtest,1)} = 0.05
        options.NRFTest (1,1) double {mustBeInteger,mustBePositive} = 1000
        options.NTestLocationsFTest (1,1) double {mustBeInteger,mustBePositive} = 10000
        options.A0 (1,1) double {mustBeNonnegative} = 0
        options.LowDensityThreshold (1,1) double {mustBeNonnegative} = 0
        options.bL (1,1) double {mustBePositive} = 50
        
        options.KFunctionDistances double = []
        options.BRandomLabelling (1,1) double {mustBeInteger,mustBePositive} = 1000
        options.distanceSpanRelationshipSign double = []
        options.alphaRandomLabelling (1,1) double {mustBePositive,mustBeLessThan(options.alphaRandomLabelling,1)} = 0.05
    end

    Nobs = size(Xobs,1);
    if(isempty(options.outDirectoryRandomBaseline))
        options.outDirectoryRandomBaseline = fullfile(pwd,"RandomBaseline");
    end
    if(~isfolder(options.outDirectoryRandomBaseline))
        mkdir(options.outDirectoryRandomBaseline)
    end
    if(isempty(options.outDirectorySHAPPE))
        options.outDirectorySHAPPE = fullfile(pwd,"SHAPPE");
    end
    if(~isfolder(options.outDirectorySHAPPE))
        mkdir(options.outDirectorySHAPPE)
    end
    if(isempty(options.outDirectoryRepulsiveRegions))
        options.outDirectoryRepulsiveRegions = fullfile(pwd,"RepulsiveRegions");
    end
    if(~isfolder(options.outDirectoryRepulsiveRegions))
        mkdir(options.outDirectoryRepulsiveRegions)
    end
    if(isempty(options.outDirectoryAttractiveRegions))
        options.outDirectoryAttractiveRegions = fullfile(pwd,"AttractiveRegions");
    end
    if(~isfolder(options.outDirectoryAttractiveRegions))
        mkdir(options.outDirectoryAttractiveRegions)
    end

    options.ShapeObs = lower(string(options.ShapeObs));
    options.shapeSynthetic = lower(string(options.shapeSynthetic));
    % Validate the shapes
    if(~((strcmp(options.ShapeObs,"square") || strcmp(options.ShapeObs,"circle")) && ...
            (strcmp(options.shapeSynthetic,"square") || strcmp(options.shapeSynthetic,"circle"))))
        errMsg = "Only 'square' or 'circle' is implemented as shapes for now.";
        error(errMsg)
    end
    if(options.RObs < 0 && options.AreaObs < 0)
        options.RObs = max(max(abs(Xobs(:,1))), max(abs(Xobs(:,2))));
        if(strcmpi(options.ShapeObs,"square"))
            options.AreaObs = (2.*options.RObs).^2;
        elseif(strcmpi(options.ShapeObs,"circle"))
            options.AreaObs = pi.*(options.RObs.^2);
        end
    elseif(options.RObs < 0)
        if(strcmpi(options.ShapeObs,"square"))
            options.RObs = sqrt(options.AreaObs)./2;
        elseif(strcmpi(options.ShapeObs,"circle"))
            options.RObs = sqrt(options.AreaObs./pi);
        end
    elseif(options.AreaObs < 0)
        if(strcmpi(options.ShapeObs,"square"))
            options.AreaObs = (2.*options.RObs).^2;
        elseif(strcmpi(options.ShapeObs,"circle"))
            options.AreaObs = pi.*(options.RObs.^2);
        end
    end

    if(isempty(options.KFunctionDistances))
        options.KFunctionDistances = linspace(5,(options.RObs)./2,50);
    end
    if(isempty(options.distanceSpanRelationshipSign))
        options.distanceSpanRelationshipSign = options.KFunctionDistances([1 end]);
    end

    if(isempty(options.types))
        options.types = unique(Xobs(:,end))';
    else
        % Check if types contains all types in XObs
        if(~all(ismember(Xobs(:,end), options.types)))
            errMsg = "All types in XObs must be present in options.types.";
            error(errMsg)
        end
    end
    options.types = options.types(:).';
    G = numel(options.types);
    NgObs = zeros(1,G);
    for g=1:G
        NgObs(1,g) = sum(Xobs(:,end)==options.types(g));
    end
    type_proportions = NgObs./Nobs;
    validTypes = NgObs>=options.Nmin;

    if(isempty(options.types_name))
        options.types_name = "type" + string(options.types);
    end

    if(ischar(options.typeOrdering) || isstring(options.typeOrdering))
        options.typeOrdering = lower(string(options.typeOrdering));
    end
    if(isstring(options.typeOrdering) && (strcmp(options.typeOrdering,"descend") || strcmp(options.typeOrdering,"ascend")))
        typesToRank = 1:G;
        [~,idx] = sort(type_proportions(typesToRank),options.typeOrdering);
        typeOrder = options.types(idx);
    elseif(isstring(options.typeOrdering))
        errMsg = "Only 'descend' or 'ascend' are implemented as type ordering function. "+...
            "Please pass directly the type ordering as a row vector if you want to use another type ordering function.";
        error(errMsg)
    else
        % Check that types in options.typeOrdering are not the same as in options.types
        typeOrder = options.typeOrdering;
        typeOrder = typeOrder(:).'; % Transform to a row vector if it is not already
        if(~isequal(sort(typeOrder),sort(options.types)))
            errMsg = "options.typeOrdering should be a permutation of options.types.";
            error(errMsg)
        end
    end

    if(ischar(options.C) || isstring(options.C))
        options.C = lower(string(options.C));
    end
    if(strcmp(options.C,"gaussian"))
        C_integerEmbedded = 1;
    elseif(strcmp(options.C,"t"))
        C_integerEmbedded = 2;
    else
        errMsg = "Only copula family 'Gaussian' or 't' or implemented for now.";
        error(errMsg)
    end

    % For reproducibility
    base_seed = 123456;
    rng(base_seed, 'twister');

    % Create a structure for the outputs
    outputs = struct();
    %%%% Step 0.1: Homogeneity Identification
    %%%%%% Step 1: Ripley's F-function    
    timerHomogeneityIdentification = tic;
    [step1Result,~] = RipleyFTest(Xobs(:,1:2),options.RObs,options.ShapeObs,'alpha',options.alphaFtest,'plotEnvelope',false,...
        'NR',options.NRFTest,'NTestLocations',options.NTestLocationsFTest,'exclusionDistance',options.exclusionDistance);
    
    LowDensityRegionsAnalysis = struct();
    if(step1Result)
        %%%% Step 2: Minimal area threshold
        thresholdInformation = struct();
        thresholdInformation.DensityThreshold = options.LowDensityThreshold;
        thresholdInformation.minArea = options.A0;
        LowDensityRegionsAnalysis = densityGridAnalysis(Xobs,options.RObs,options.bL,options.ShapeObs,thresholdInformation);
        %%% How many low density regions sections?
        N_LowDensity = LowDensityRegionsAnalysis.num_patches;
        inhomogeneous = (N_LowDensity>0);
    else
        inhomogeneous = step1Result;
    end
    
    elapsedTimeHomogeneityIdentification = toc(timerHomogeneityIdentification);

    outputs.HomogeneityIdentificationElapsedTime = elapsedTimeHomogeneityIdentification;
    outputs.HomogeneityStatus = ~inhomogeneous;
    outputs.LowDensityRegionsAnalysis = LowDensityRegionsAnalysis;

    %%%% Step 0.2: Evaluate clustering patterns using random labelling
    Ndistances = numel(options.KFunctionDistances);
    HvalsMatrix = nan(G,G,Ndistances);
    weightMatrices = nan(Nobs,Nobs,Ndistances);
    timerRandomLabelling = tic;
    %%% 1. Evaluate the K-function
    for gg1=1:G
        g1 = options.types(gg1);
        if(NgObs(1,gg1)>=options.Nmin)
            for gg2=1:G
                g2 = options.types(gg2);
                if(NgObs(1,gg2)>=options.Nmin)
                    [Hvals,marks,pairwise_weights] = ripleyK(Xobs, g1, g2, options.RObs, options.KFunctionDistances,'H',options.ShapeObs);
                    HvalsMatrix(gg1,gg2,:) = Hvals;
                    % Fill weights matrix
                    for d=1:Ndistances
                        weightMatrices(Xobs(:,end)==g1,Xobs(:,end)==g2,d) = pairwise_weights(d).wij;
                    end
                end
            end
        end
    end
    
    %%% 2. Fill the weight matrix for dropped types
    for gg1=1:G
        g1 = options.types(gg1);
        for gg2=1:G
            g2 = options.types(gg2);
            if((NgObs(1,gg1)<options.Nmin || NgObs(1,gg2)<options.Nmin) && NgObs(1,gg1)>0 && NgObs(1,gg2)>0)
                [~,~,pairwise_weights] = ripleyK(Xobs, g1, g2, options.RObs, options.KFunctionDistances,'H',options.ShapeObs);
                % Fill weights matrix
                for d=1:Ndistances
                    weightMatrices(Xobs(:,end)==g1,Xobs(:,end)==g2,d) = pairwise_weights(d).wij;
                end
            end
        end
    end
    
    %%% 3. Generate B random permutations of types and evaluate K-function for
    %%% each permutation
    nullDistributionsH = nan(G, G, Ndistances, options.BRandomLabelling);
    
    lambdaI = NgObs'./options.AreaObs; % Gx1 matrix
    lambdaJ = repmat(NgObs./options.AreaObs,G,1); % GxG matrix
    for g=1:G
        lambdaJ(g,g) = (NgObs(1,g)-1)./options.AreaObs; % Univariate case
    end
    normMatrix = (lambdaI.*lambdaJ).*options.AreaObs; % GxG matrix
    
    for b = 1:options.BRandomLabelling
        randomIndices = randperm(Nobs);
        newmarks = Xobs(randomIndices, end);
        M = (newmarks == options.types); % NobsxG matrix where element (i,j) is 1 if point i is type j and 0 otherwise
        for rr = 1:Ndistances
            S = (M.'*(weightMatrices(:,:,rr)*M))./normMatrix;
            nullDistributionsH(:,:,rr,b) = sqrt(S./pi)-options.KFunctionDistances(rr);
        end
    end
    
    % Apply Nmin
    nullDistributionsH(~validTypes,:,:,:) = NaN;
    nullDistributionsH(:,~validTypes,:,:) = NaN;
    
    [~,idx1] = min(abs(options.KFunctionDistances-options.distanceSpanRelationshipSign(1)));
    [~,idx2] = min(abs(options.KFunctionDistances-options.distanceSpanRelationshipSign(2)));
    randomEnvelopeBounds = nan([G G 3 numel(options.KFunctionDistances(idx1:idx2))]);
    
    corrMatrix = nan.*ones(G,G);
    
    for gg1=1:G
        g1 = options.types(gg1);
        if(NgObs(1,gg1)>=options.Nmin)
            for gg2=1:G
                g2 = options.types(gg2);
                if(NgObs(1,gg2)>=options.Nmin)
                    nullDistributionValues = reshape(nullDistributionsH(gg1,gg2,:,:),[Ndistances options.BRandomLabelling])';
                    mu = mean(nullDistributionValues, 1,'omitnan');                 
                    sd = std(nullDistributionValues, 0, 1,'omitnan');              
                    sd(sd==0) = eps;
                    RipleyFunction = reshape(HvalsMatrix(gg1,gg2,:),[1,Ndistances]);
                    Zobs = (RipleyFunction - mu) ./ sd;
                    Zsim = (nullDistributionValues - mu) ./ sd;
                    
                    % Two-sided global test (max absolute deviation)
                    Tobs = max(abs(Zobs));           
                    Tsim = max(abs(Zsim), [], 2);        
                    
                    pval = (sum(Tsim >= Tobs)+1)./(options.BRandomLabelling + 1);
                    corrMatrix(gg1,gg2) = pval;
    
                    % Lower bounds
                    randomEnvelopeBounds(gg1,gg2,1,:) = min(nullDistributionValues(:,idx1:idx2));
                    % Upper bounds
                    randomEnvelopeBounds(gg1,gg2,2,:) = max(nullDistributionValues(:,idx1:idx2));
                    % Function
                    randomEnvelopeBounds(gg1,gg2,3,:) = RipleyFunction(idx1:idx2);
                end
            end
        end
    end
    
    % Benjamini-Hochberg Correction
    idx_not_nan = find(~isnan(corrMatrix));
    pvals_flat = corrMatrix(idx_not_nan);
    [adjusted_p_values] = mafdr(pvals_flat,'BHFDR',true);
    corrMatrix(idx_not_nan) = adjusted_p_values;
    
    % Are you significantly above or under?
    aboveMatrix = zeros(G,G);
    belowMatrix = zeros(G,G);
    for gg1=1:G
        g1 = options.types(gg1);
        if(NgObs(1,gg1)>=options.Nmin)
            for gg2=1:G
                g2 = options.types(gg2);
                if(NgObs(1,gg2)>=options.Nmin)
                    if(corrMatrix(gg1,gg2)<=options.alphaRandomLabelling)
                        LBs = reshape(randomEnvelopeBounds(gg1,gg2,1,:),[1,numel(options.KFunctionDistances(idx1:idx2))]);
                        UBs = reshape(randomEnvelopeBounds(gg1,gg2,2,:),[1,numel(options.KFunctionDistances(idx1:idx2))]);
                        K = reshape(randomEnvelopeBounds(gg1,gg2,3,:),[1,numel(options.KFunctionDistances(idx1:idx2))]);
                        Aabove = sum(max(K-UBs,0));
                        Abelow = sum(max(LBs-K,0));
                        if(Aabove>Abelow)
                            aboveMatrix(gg1,gg2) = aboveMatrix(gg1,gg2)+1;
                        elseif(Aabove<Abelow)
                            belowMatrix(gg1,gg2) = belowMatrix(gg1,gg2)+1;
                        end % end of if(Aabove>Abelow) ... elseif ... 
                    end % end of if(corrMatrix(gg1,gg2)<=alpha)
                end % end of if(Ng(1,gg2)>=Nmin)
            end % end of for gg2=1:G
        end % end of if(Ng(1,gg1)>=Nmin)
    end % End of for gg1=1:G
    
    elapsedTimeRandomLabelling = toc(timerRandomLabelling);

    outputs.RandomLabellingTestElapsedTime = elapsedTimeRandomLabelling;
    outputs.HvalsMatrix = HvalsMatrix;
    outputs.weightMatrices = weightMatrices;
    outputs.nullDistributionsH = nullDistributionsH;
    outputs.aboveMatrix = aboveMatrix;
    outputs.belowMatrix = belowMatrix;
    outputs.corrMatrix = corrMatrix;

    %%%% SHAPPE Reproductions
    timerSHAPPE = tic;
    
    % Get distance distributions from observed data
    [distributions] = getDistanceDistributions(Xobs,options.types,options.types_name,options.L); 
    convergenceOutputs_allReplicates = zeros(options.Nreps,4,G);
    
    % Loop through replicates
    parfor(nr=1:options.Nreps,options.parallelLoops)
        rng(base_seed + 1000*nr, 'twister');
    
        % Step 0.3: Generate and save random baseline
        randomBaselineFilename = fullfile(options.outDirectoryRandomBaseline,"rep"+string(nr)+".csv");
        attractiveRegionsFilename = fullfile(options.outDirectoryAttractiveRegions,"rep"+string(nr)+".csv");
        repulsiveRegionsFilename = fullfile(options.outDirectoryRepulsiveRegions,"rep"+string(nr)+".csv");
        [randomBaselineConfiguration,R,repulsiveRegions,attractiveRegions,Ng] = ...
            ECSRGenerator(options.N,type_proportions,options.N./Nobs,options.RObs,options.shapeSynthetic,...
            options.ShapeObs,options.exclusionDistance,options.types,~inhomogeneous,options.attract,LowDensityRegionsAnalysis);
        writematrix(randomBaselineConfiguration,randomBaselineFilename);
        writematrix(attractiveRegions,attractiveRegionsFilename);
        writematrix(repulsiveRegions,repulsiveRegionsFilename);
    
        if(options.attract && options.useRegions)
            Regions = attractiveRegions;
        elseif(options.useRegions)
            Regions = repulsiveRegions;
        else
            Regions = [];
        end
    
        % Create array for convergence information for this replicate
        convergenceOutput = zeros(4,G);
    
        % Get distance distributions from ECSR
        [distributionsRandom] = getDistanceDistributions(randomBaselineConfiguration,options.types,options.types_name,options.L);
    
        % Step 1: Start from ECSR configuration
        SHAPPEconfig = randomBaselineConfiguration;
    
        % Information of types in order they will be placed
        typedata = struct();
        for g=1:G
            typedata(g).ObsDist = distributions.(options.types_name(options.types==typeOrder(g)));
            typedata(g).RandomDist = distributionsRandom.(options.types_name(options.types==typeOrder(g)));
        end % End of for g=1:G
    
        samplingindexes = zeros(2,G);
        for gg_idx=1:G
            samplingindexes(1,gg_idx) = options.L*(gg_idx-1)+1;
            samplingindexes(2,gg_idx) = (options.L*(gg_idx));
        end
    
        to_exclude = ~validTypes;
        number_points_placed = 0;
        placedTypes = [];
        macroIdxSampled = [];
    
        for gg=1:G
            current_type = typeOrder(gg);
            current_type_idx = find(options.types==current_type);
            % Find  points
            points_index = find(SHAPPEconfig(:,end)==current_type);
            Ncurrent = Ng(current_type_idx); % total number of cells of current type in the simulation
            
            % Step 2: Minimal threshold to perform iterations
            if(~to_exclude(current_type_idx))
                dist = typedata(gg).ObsDist;
                distRandom = typedata(gg).RandomDist;
    
                % Place points randomly but make sure they are
                % at least h units away from other points
                SHAPPEconfig = placeWithExclusion(SHAPPEconfig, points_index, placedTypes, Regions, options, R,1000);

                % Step 3: Is the type spatially correlation to at least one
                % type amongst self-type and already placed types?
                sampling_index= (samplingindexes(1,current_type_idx):samplingindexes(2,current_type_idx));
                if(~isempty(placedTypes))
                    for pct=1:numel(placedTypes)
                        pct_idx = find(options.types==placedTypes(pct));
                        sampling_index = [sampling_index,(samplingindexes(1,pct_idx):samplingindexes(2,pct_idx))];
                    end % end of for pct=1:numel(placedTypes)
                end % end of if(~isempty(placedTypes))
                sampling_index = sort(sampling_index);
    
                uncorrelated_indexes = [];
                to_exclude_current = find(to_exclude);
                % Check autocorrelation status
                if(corrMatrix(current_type_idx,current_type_idx)>options.SHAPPEalpha || isnan(corrMatrix(current_type_idx,current_type_idx)))
                    to_exclude_current = unique([to_exclude_current,current_type_idx]);
                    sampling_index = setdiff(sampling_index,(samplingindexes(1,current_type_idx):samplingindexes(2,current_type_idx)));
                    uncorrelated_indexes = [uncorrelated_indexes,(samplingindexes(1,current_type_idx):samplingindexes(2,current_type_idx))];
                end
                % Check correlation status to other cell types
                if(~isempty(placedTypes))
                    for pct=1:numel(placedTypes)
                        pct_idx = find(options.types==placedTypes(pct));
                        if(corrMatrix(current_type_idx,pct_idx)>options.SHAPPEalpha || ...
                                isnan(corrMatrix(current_type_idx,pct_idx)) || ...
                                ~validTypes(pct_idx))
                            to_exclude_current = unique([to_exclude_current,pct_idx]);
                            sampling_index = setdiff(sampling_index,(samplingindexes(1,pct_idx):samplingindexes(2,pct_idx)));
                            uncorrelated_indexes = [uncorrelated_indexes,(samplingindexes(1,pct_idx):samplingindexes(2,pct_idx))];
                        end
                    end
                end
    
                compare_with = [placedTypes,current_type];
                comparisonWeights = ones(1,G);
                if(~isempty(to_exclude_current))
                    for ee=1:numel(to_exclude_current)
                        comparisonWeights(to_exclude_current(ee)) = options.wNoCorr;
                    end
                end
    
                if(isempty(sampling_index))
                    % STEP 5: Impose environmental constraints
                    [SHAPPEconfig,macroIdxSampled] = EnvironmentalConstraints(SHAPPEconfig,Xobs,options.bLEnv,...
                        options.RObs,R,current_type,[placedTypes,current_type],...
                        options.exclusionDistance,options.ShapeObs,options.shapeSynthetic,macroIdxSampled);
                    achievedEndCriterion = 4;
                else
                    % STEP 4
                    %%% Sample distances with a copula
                    ObsMatrix = [];
                    for g=1:G
                        ObsMatrix = [ObsMatrix;dist.(options.types_name(g)).data];
                    end
                    ObsMatrix = ObsMatrix';
                    O = ObsMatrix(:, sampling_index);
                    [Rtau,~] = corr(O, 'type', 'Kendall', 'Rows', 'pairwise');
                    % Map Kendall's Tau to latent Gaussian correlation
                    rho = sin( (pi/2) * Rtau );
                    % Avoid borderline cases
                    rho(isnan(rho)) = 0; %XX
                    rho = max(min(rho, 0.999), -0.999);
                    % Ensure symmetry - may have been lost because of
                    % floating points
                    rho = (rho+rho')./2;
                    % Make all elements on the diagonal be equal to one
                    rho = rho - diag(rho).*eye(size(rho,1)) + eye(size(rho,1));
                    % Nearest correlation matrix in case rho is not PD
                    % The if loop is in case nearcorr does not work
                    if(exist('nearcorr','file')==2)
                        rho = nearcorr(rho);
                        rho = (rho + rho.')/2;
                        rho(1:size(rho,1)+1:end) = 1;
                    else
                        [V,D] = eig((rho+rho.')/2);
                        d = diag(D);
                        eps_bump = 1e-10;
                        if min(d) < eps_bump
                            D = diag(max(d, eps_bump));
                            rho = V*D*V.';
                            rho = (rho + rho.')/2;
                            rho(1:size(rho,1)+1:end) = 1;
                        end % End of if min(d) < eps_bump
                    end % End of if exist('nearcorr','file') == 2 ... else ... end
    
                
                    if(length(sampling_index)>1)
                        if(C_integerEmbedded>1)
                            % Fit the degree of freedom nu in a t copula
                            %%% Transform matrix of observations O to uniform
                            %%% marginals
                            O = ObsMatrix(:,sampling_index);
                            unif_marg = tiedrank(O)./(size(O,1)+1);
                            nu_vals = 2:0.5:50;                 
                            logL = zeros(size(nu_vals));
                            for i = 1:numel(nu_vals)
                                logL(i) = sum(log(copulapdf('t', unif_marg, rho, nu_vals(i))));
                            end
                            [~, idx] = max(logL);
                            nu = nu_vals(idx);
                        end
                    end
                    target_distances = zeros(options.N,G*options.L);
                    
                    copulaInformation = struct();
                    copulaInformation.rho = rho;
        
                    if(length(sampling_index)>1)
                        copulaInformation.kernel = options.C;
                        if(C_integerEmbedded>1)
                            copulaInformation.nu = nu;
                        end
                    else
                        copulaInformation.kernel = "None";
                    end
        
                    for nC=(points_index)'
                        target_distances(nC,:) = sampleTargetDistance(options.L,dist,distRandom,copulaInformation,options.types,...
                            samplingindexes,options.types_name,sampling_index,uncorrelated_indexes,options.exclusionDistance);
                    end
        
                    % Perform iterative point placement
                    if(~options.useRegions)
                        regionArgs = {};
                    elseif(options.attract)
                        regionArgs = {'attractiveRegions',Regions};
                    else
                        regionArgs = {'repulsiveRegions',Regions};
                    end
                    [SHAPPEconfig,~,achievedEndCriterion,~,~,macroIdxSampled] = IterativePointPlacement(SHAPPEconfig,Xobs,target_distances,...
                        options.types,current_type,compare_with,comparisonWeights,placedTypes,options.RObs,R,options.ShapeObs,options.shapeSynthetic,...
                        'maxIterations',options.MaxN,'phi',options.phi,'phiMinimum',options.phiMinimum,'rW',options.rW,'L',options.L,...
                        'a',options.a,'lambda',options.lambda,'macroIdxSampled',macroIdxSampled,'macroBinLength',options.bLEnv,...
                        'ghostDistance',options.ghostDistance,'exclusionDistance',options.exclusionDistance,'Nenv',options.nEnv,...
                        'dMed',options.dMed,'d90',options.d90,'dNet',options.dNet,'nEq',options.nEq,regionArgs{:});
                end
                convergenceOutput(achievedEndCriterion,current_type_idx) = 1;
                number_points_placed = number_points_placed+Ncurrent;
                placedTypes = [placedTypes,current_type];
            else
                if(~isempty(points_index))
                    SHAPPEconfig = placeWithExclusion(SHAPPEconfig, points_index, placedTypes, Regions, options, R,1000);
                    % Impose environmental contraints
                    [SHAPPEconfig,macroIdxSampled] = EnvironmentalConstraints(SHAPPEconfig,Xobs,options.bLEnv,...
                        options.RObs,R,current_type,[placedTypes,current_type],...
                        options.exclusionDistance,options.ShapeObs,options.shapeSynthetic,macroIdxSampled);
                end
                number_points_placed = number_points_placed+Ncurrent;
                placedTypes = [placedTypes,current_type];
                convergenceOutput(:,current_type_idx) = -1.*ones(4,1);
            end % End of if(any(to_exclude==current_cell_type))
    
        end % End of for gg=1:G
        convergenceOutputs_allReplicates(nr,:,:) = convergenceOutput;
        writematrix(SHAPPEconfig,fullfile(options.outDirectorySHAPPE,"rep"+string(nr)+".csv"))
    end % end of parfor nr=1:Nreps
    elapsedTimeSHAPPE = toc(timerSHAPPE);
    outputs.elapsedTimeSHAPPE = elapsedTimeSHAPPE;
    outputs.convergenceOutputs_allReplicates = convergenceOutputs_allReplicates;

end