function [densityMaps,step,pointsIndexes] = buildDensityMap(conf,R,types,varargin)
    
    p = inputParser;    
    addRequired(p,'conf');
    addRequired(p,'R');
    addRequired(p,'types');
    addParameter(p,'nbin',[]);
    addParameter(p,'step',[]);
    addParameter(p,'units',"probability");
    
    parse(p, conf, R, types, varargin{:});

    nbin = p.Results.nbin;
    step = p.Results.step;
    units = p.Results.units;

    if(isempty(step) && isempty(nbin))
        nbin = 25;
        step = 2.*R./nbin;
    elseif(isempty(step))
        step = 2.*R./nbin;
    elseif(isempty(nbin))
        nbin = round(2.*R./step);
        step = 2.*R./nbin;
    else
        msg = "You must enter either no value for nbin and step, or a value for one of them, but not both.";
        error(msg)
    end

    % We wand positions between [0,2*R]
    if(any(conf(:,1)<0) || any(conf(:,2)<0))
        conf(:,1) = R+conf(:,1);
        conf(:,2) = R+conf(:,2);
    end

    densityMaps = zeros(numel(types),nbin,nbin);
    pointsIndexes = struct;
    pointsIndexes(numel(types)) = struct;
    for t=types(:).'
        idxTypes = find(conf(:,end)==t);
        pointsIndexes(types==t).bins(nbin.*nbin).agents = [];
        for nc=idxTypes(:).'
            x = conf(nc,1);
            y = conf(nc,2);
            xbin = min(nbin,max(ceil(x/step),1));
            ybin = min(nbin,max(ceil(y/step),1));
            densityMaps(types==t,xbin,ybin) = densityMaps(types==t,xbin,ybin)+1;
            % find the bin within which the point falls
            bins_idx = sub2ind([nbin, nbin], xbin, ybin);
            pointsIndexes(types==t).bins(bins_idx).agents = [pointsIndexes(types==t).bins(bins_idx).agents,nc];
        end % End of for nc=idxTypes
        if(strcmp(units,"probability") && ~isempty(idxTypes))
            densityMaps(types==t,:,:) = densityMaps(types==t,:,:)./numel(idxTypes);
        end
    end % End of for t=1:length(types)

end