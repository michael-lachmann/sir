%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/zwdu/Documents/GitHub/sir/pholme_sir/AustinSmall10Net10.csv
%
% Auto-generated by MATLAB on 30-Apr-2020 17:10:59
function AustinSmall10Net10 = readNetFile(filename)
%% Setup the Import Options and import the data
% # network files have 5 columns:
% #   FROM TO GR1 GR2 weight
% # which gives one contact edge. weight is an integer giving the index of the weight you specify the actual value in the run
% # GR1 and GR2 are the groups (age groups currently) of the nodes

opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["FROM", "TO", "GR1", "GR2", "weight"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
AustinSmall10Net10 = readtable(filename, opts);


%% Clear temporary variables
clear opts