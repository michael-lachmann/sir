
%% initialize all nodes to susceptible
S = ones(len_Nodes,1);
I = S*0; R = I;

%% initialize list of Recovery time
RecT_Max = 1000;
RecT = RecT_Max*S; 

%% set a random source node i into SIR model
ini_Node = randi(length(S));
S(ini_Node) = 0;
I(ini_Node) = 1;
R(ini_Node) = 1;

%% rate of S to I and I to R
rateSI = 0.15;
rateIR = 0.1;

%% label of different status
S_lab = 1;
I_lab = 2;
R_lab = 3;
Labs = [1:3];
Labs_Str = {'S', 'I', 'R'};

%% initialize the heap, ranking by time in decreasing order
% heap.queue with four cols: Node id, From status, To status, time gap
% heap.nq: number of events in heap.queue.
heap.queue = zeros(100, 4);
heap.nq = 0;



%% 
NEVER=10^4;