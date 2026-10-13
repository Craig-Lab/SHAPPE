function new_config = mapIMC_to_ABM(config)

    imc_gbm = 1;
    imc_cd8 = 2;
    imc_stroma = 3;
    imc_m1 = 4;
    imc_m2 = 5;
    imc_m0 = 6;
    
    abm_gbm = 2;
    abm_cd8 = 1;
    abm_stroma = 3;
    abm_m1 = 5;
    abm_m2 = 6;
    abm_m0 = 4;
    
    mapping = [imc_gbm,abm_gbm;...
        imc_cd8,abm_cd8;...
        imc_stroma,abm_stroma;...
        imc_m1,abm_m1;...
        imc_m2,abm_m2;...
        imc_m0,abm_m0];
    
    new_config = config;
    for c=1:size(mapping,1)
        c_new = mapping(c,1);
        c_index = find(config(:,end)==c_new);
        if(~isempty(c_index))
            new_config(c_index,end) = mapping(c,2);
        end
    end
end