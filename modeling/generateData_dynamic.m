function [data] = generateData_dynamic(numChoices, numAgents, numAtts, numValuesPerAtt, ...
    params, param_structure, param_structure_step, sign_optionvals)

numChoicesPerReset = 10;

lbs1 = [param(1).lb repelem(param(2).lb, numAtts) repelem(param(numAtts+2).lb, numAtts)]; % full weights
ubs1 = [param(1).ub repelem(param(2).ub, numAtts) repelem(param(numAtts+2).ub, numAtts)];
lbs2 = [param(1).lb repelem(param(2).lb2, numAtts) repelem(param(numAtts+2).lb, numAtts)]; % signed weights
ubs2 = [param(1).ub repelem(param(2).ub2, numAtts) repelem(param(numAtts+2).ub, numAtts)];

for agent = 1:numAgents

    % generate first set of options
    options = signOptions(randi([1 numValuesPerAtt], numAtts, 2, numChoicesPerReset) / numValuesPerAtt, ...
        sign_optionvals, false);
    avail_atts = generateAvailAtts(numChoicesPerReset, numAtts, 1, 7);

    cur_choice_index = 1:numChoicesPerReset;

    data(agent).N = numChoices;
    data(agent).params = params(agent,:);

    data(agent).options = options;
    data(agent).avail_atts(cur_choice_index,:) = avail_atts;
    [data(agent).choices(cur_choice_index)] = makeChoice(data(agent).params, options, avail_atts);

    for i = 2:(numChoices/numChoicesPerReset)
        % do fitting
        % WAD
        WAD_post = @(x) -post_WAD(x,data(agent),false,param_structure,numAtts);
        [results_WAD.x,~] = ga(WAD_post, length(param_structure),[],[],[],[],lbs1,ubs1,[]);
        results_WAD.likfun = @(x,d) getLogLik_step(x, d, false);

        % WP
        WP_post = @(x) -post_WAD(x,data(agent),true,param_structure,numAtts);
        [results_WP.x,~] = ga(WAD_post, length(param_structure),[],[],[],[],lbs1,ubs1,[]);
        results_WP.likfun = @(x,d) getLogLik_step(x, d, true);

        % EW
        EW_post = @(x) -post_WAD(x,data(agent),true,param_structure,numAtts);
        [results_WP.x,~] = ga(WAD_post, length(param_structure),[],[],[],[],lbs1,ubs1,[],(numParams-numAtts+1):numParams);
        results_WP.likfun = @(x,d) getLogLik_step(x, d, true);

        % TAL
        TA:_post = @(x) -post_WAD(x,data(agent),true,param_structure,numAtts);
        [results_WP.x,~] = ga(WAD_post, length(param_structure),[],[],[],[],lbs1,ubs1,[],(numParams-numAtts+1):numParams);
        results_WP.likfun = @(x,d) getLogLik_step(x, d, true);

        % WAD_step
        WAD_post_step = @(x) -post_WAD_step(x,data(agent),false,param_structure_step,numAtts);
        [results_WAD_step.x,~] = ga(WAD_post_step, length(param_structure_step),[],[],[],[],lbs1,ubs1,[],(numParams-numAtts+1):numParams);
        results_WAD_step.likfun = @(x,d) getLogLik_step(x, d, false);


        % generate options
        optimizer = @(opts);
    
        % simulate choices
        cur_choice_index = (1+numChoicesPerReset*(i-1)):(numChoicesPerReset*i);
    
        for agent = 1:numAgents
            data(agent).N = numChoicesPerReset*i;
            data(agent).options(:,:,cur_choice_index) = options(:,:,:,agent);
            data(agent).avail_atts(cur_choice_index,:) = avail_atts(:,:,agent);
            [data(agent).choices(cur_choice_index)] = makeChoice(data(agent).params, data(agent).options, data(agent).avail_atts);
        end
    end
end