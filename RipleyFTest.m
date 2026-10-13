function [h,pval,varargout] = RipleyFTest(pointPattern,R,shape,varargin)
    p = inputParser;    
    addRequired(p,'pointPattern');
    addRequired(p,'R');
    addRequired(p,'shape');
    addParameter(p,'alpha',0.05);
    addParameter(p,'plotEnvelope',false);
    addParameter(p,'NR',1000);
    addParameter(p,'NTestLocations',10000);
    addParameter(p,'returnInformation',false);
    addParameter(p,'exclusionDistance',0);
    
    parse(p, pointPattern, R,shape, varargin{:});

    alpha = p.Results.alpha;
    plotEnvelope = p.Results.plotEnvelope;
    NR = p.Results.NR;
    NTestLocations = p.Results.NTestLocations;
    returnInformation = p.Results.returnInformation;
    exclusionDistance = p.Results.exclusionDistance;

    if(strcmp(shape,"circle"))
        dist_max = 2.*R;
    elseif(strcmp(shape,"square"))
        dist_max = sqrt(2.*(2.*R).^2);
    else
        msg = "Shape must be circle or square.";
        error(msg)
    end % End of if(strcmp(shape,"circle")) ... elseif(strcmp(shape,"square")) ... else ... end
    distances = 5:5:dist_max;
    F = ripleyF(pointPattern,R,shape,distances,'NTestLocations',NTestLocations);
    O = distances(find(F>=1,1,'first'));

    % Generate null distribution
    nullDistributionMax = zeros(1,NR);
    nullDistribution = zeros(NR,length(distances));
    for r=1:NR
        randomConfig = ECSRGenerator(size(pointPattern,1), 1, 1, R, shape, shape, exclusionDistance, 1, true);
        Fnull = ripleyF(randomConfig(:,1:2),R,shape,distances,'NTestLocations',NTestLocations);
        nullDistributionMax(1,r) = distances(find(Fnull>=1,1,'first'));
        nullDistribution(r,:) = Fnull(:).';
    end % End of for r=1:NR

    rmax_idx = find(distances<=max(max(nullDistributionMax(1,:)),O),1,'last');
    if(isempty(rmax_idx))
        rmax_idx = length(distances);
    end % End of if(isempty(rmax_idx))
    nullDistribution = nullDistribution(:,1:rmax_idx);
    F = F(1:rmax_idx);
    distances = distances(1:rmax_idx);

    % Calculate p value and test result
    mu = mean(nullDistribution,1);
    sd = std(nullDistribution,1);
    % Avoid blow ups because of very small deviation
    keep = sd > 1e-6;
    F = F(keep);
    nullDistribution = nullDistribution(:,keep);
    mu = mu(keep);
    sd = sd(keep);
    distances = distances(keep);
    Zobs = (mu - F) ./ sd;
    Zsim = (mu - nullDistribution) ./ sd;
    Tobs = max(Zobs);
    Tsim = max(Zsim, [], 2); 
    pval = (sum(Tsim >= Tobs) + 1) / (NR + 1);
    h = (pval<=alpha);

    argOutNumber = 1;
    if(plotEnvelope)
        fig = figure;
        hold on
        plot_areaerrorbar(nullDistribution(:,:),[],distances)
        plot(distances,F,'LineWidth',1.5,'Color','red')
        varargout{argOutNumber} = fig;
        argOutNumber = argOutNumber+1;
    end
    if(returnInformation)
        info = struct();
        info.nullDistribution = nullDistribution;
        info.F = F;
        info.distances = distances;
        varargout{argOutNumber} = info;
        argOutNumber = argOutNumber+1;
    end
end