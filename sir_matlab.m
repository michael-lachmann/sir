%% // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  // code for SIR on networks by Petter Holme (2018)
%
%% Michael's comments when zw checked on Apr 28, 2020:
% // One question is if to go through the whole disease progression when we pick an individual up.
% // Is the heap just for infections or for any events?
% // Events other than infections are independent. No one cares how it progresses.
% // Therefore, putting them on a heap is a waste. You don't care if they occured before or after someone else, they do not cause any further events in other individuals.

% // However, in the Meyer lab mode, you are still infectious when you do transitions. So mayeb then it does make sense to put them on the heap?
% // So maybe the rule should be do as little ass possible so you don't have to do things twice?

% // If we knew when the next change in the network would be, we could compare to it. But is it worth it to compare all the time, when it is so rare?
g.w_n = 0 ;
g.av_deg = 0.0 ;
g.ngroup = 1 ;
NAVG=10; % // number of runs for averages
NSTATES =(256);
NQUART = 65536; %0x10000;
% RANDOM_QUART = &rand;
EPSILON = (1e-10);
NEVER = (10^6);
MAXRAND = (4294967295.0);

UINT_MAX = 4294967295;
I_OR_R = (UINT_MAX - 1);
NONE = UINT_MAX;

set_parameters;

% // basically we get to a node, and change its state, make its time earlier
% // or later.
% // in special cases, we infect other nodes, and make their time earlier or later
% // (later only in the special case of changing parameters)
% // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% // this routine runs one SIR outbreak from a random seed node
%  init_mats_seir_full(...) init_mats_seir_full( &(struct named_seir_args){__VA_ARGS__})
% // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - 
%%     sir f,AustinSmall10Net10.csv R,63 S,B0.03 w,20,1,1,1 w,50,1,1,1 w,200,0.25,1,0
% 	case 'f': % // read network with groups from file
% 		fp = fopen(argv(k)+2, "r"); % // skip over 'F,'
% 		read_data5(fp);
g.state = 100; % random seed
g.beta = 0.03; %
pos = 2 ;
r = 0 ; npos = 0 ;

% //init_mats_seir( g.beta, g.gamma, g.sigma) ;
args.del=g.ngroup;
args.beta = g.beta; 
args.gammaA = g.gamma_a; 
args.gammaH= g.gamma_h; 
args.gammaY = g.gamma_y; 
args.sigma = g.sigma;
args.omegaY = g.omega_y; 
args.omegaA = g.omega_a(i); 
args.omegaE = g.omega_e(i); 
args.omegaH = g.omega_h ;
args.sigma = g.sigma; 
args.mu= g.mu; 
args.nu= g.nu(i); 
args.pi = g.pi(i); 
args.tau = g.tau; 
args.eta= g.eta;
for( i=1:g.ngroup)
    args.S=i;
    g = init_mats_seir_full(args, NQUART, g, MAXRAND, NEVER);
end


g.w_time = 100;
g.w_val = ones(100,1);
g.w_n = 20;
g.nweight = 10;
g.weight = ones(100,1);

for T_i=1:g.w_n
    for( i=2:g.nweight)
        g.weight(i) = (g.w_val(  T_i*(g.nweight) +i ) ) ;
    end
    if( T_i==0)
        nS = 10;
        g = seir_init(nS, NSTATES, g);
    else
        reinfect_s_I_neighbours() ;
    end
    
    last_now = g.now ;
    while( g.now < g.w_time(T_i) )
        epi_timestep;
        if (g.now - last_now > 1)
            last_now = g.now ;
        end
    end
end
% // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
