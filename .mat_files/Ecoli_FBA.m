%FBA simulation code
%load the model from the .mat file on the BiGG database website
model = load(('iJO1366.mat')).iJO1366;
%initialize vectors​
effects = zeros(length(model.rxns), 1);
ratios = zeros(length(model.rxns), 1);
ratekos = zeros(length(model.rxns), 1);
ratewts = zeros(length(model.rxns), 1);

for i = 1:length(model.rxns)
    [ratio, rateKO, rateWT, hasEffect, deleted_rxn, flux_solution] = singleRxnDeletion(model, 'FBA', model.rxns(i));
    effects(i) = hasEffect;
    ratios(i) = ratio;
    ratekos(i) = rateKO;
    ratewts(i) = rateWT;
    rxn(i) = deleted_rxn;
    %sprintf("%s %f", int2str(hasEffect), ratio)
end
rxn = rxn'; %transpose​ rxn into a column vector
had_effect = find(effects);
num_effect = length(had_effect);

%see which nodes have an effect greater than 95%
rate_ratio = zeros(length(ratekos), 1);
for i = 1:length(ratekos)
    rate_ratio(i) = ratekos(i)/ratewts(i);
end

%find which indices have an effect greater than 95%
effect95 = zeros(length(rate_ratio), 1);
for i = 1:length(rate_ratio)
    if rate_ratio(i)<0.95
        effect95(i) = 1;
    else effect95(i) = 0;
    end
end
indices = find(effect95);

%combine into a single dataframe
essential95 = [rxn, num2cell(effect95)]; %first col is rxn name, second col is whether it is essential or not
writecell(essential95, "EssentialRxns.csv"); %convert to csv

%identify the reactions that had a 0.95 effect by name from reaction list
for i = 1:length(indices)
    essentialRxns(i) = rxn(indices(i));
end
essentialRxns=essentialRxns'; %transpose so it is a column vector
%make sure to index using essentialRxns using char() in future

%Now since we don't know how they limited the original 2583 reactions to
%   1251 reactions, look at which of the essential reactions you just
%   found that are included in those 1251 reactions

adjMat=readtable('AdjData_ecoli.xlsx'); %read in ajacency matrix of 1251 rxns
rxn_list=adjMat(:,1); %take first column as reaction (node) names
rxn_list=table2array(rxn_list); %convert to array from table
final_essential = zeros(length(rxn_list), 1);
for i=1:length(rxn_list)
    currRxn = rxn_list(i);
    i;
    final_essential(i)= ismember(currRxn, essentialRxns);
end

%graph metrics (except for cascade number, which is in a separate matlab script)
%   were calculated as part of the group project, in python
