%load model
model = load(('iJO1366.mat')).iJO1366;

%load adjacency matrix from xlsx file
adjMat=readtable('AdjData_ecoli.xlsx');

%take first column as reaction (node) names
rxn_list=adjMat(:,1);
rxn_list=table2array(rxn_list);

%now remove that first column from the adjacency matrix
adjMat=adjMat(:,2:1252);
adjMat=table2array(adjMat); %convert table to an array of doubles

%create a vector that will later be filled by cascade numbers of each node
cascadeNums=zeros(length(rxn_list),1);

%for each rxn, get the cascade number and store it in cascadeNums
for i=1:length(rxn_list)
    affected_nodes=cascade(adjMat,i);
    cascadeNums(i,1)=length(affected_nodes); %length of vector containing all nodes that are affected
end

%it's working!!!

%Make sure to run Ecoli_FBA script before running this one
%Then, concatenate the final variables into two dataframes, and export
%   them into R for merging, plotting, and statistics

final_cascade=[rxn_list, num2cell(cascadeNums)];
writecell(final_cascade,'CascadeNums.csv');