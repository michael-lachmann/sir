
function g = init_mats_seir_full(args, NQUART, g, MAXRAND, NEVER)
% 	WHY?
E = args.S+  args.del,
IA= args.S+2*args.del,
IY= args.S+3*args.del,
IH= args.S+4*args.del,
R = args.S+5*args.del,
D = args.S+6*args.del,
S = args.S;

% // S
for(i=1:NQUART)
    g.time_dist (S,i) = myQexp( (1.0+i)/NQUART ) / args.beta ;
end

g.state_trans_p(S,1) = MAXRAND ;  g.state_trans(S,1) = E ;

g.state_infect_rate(S) = 0.0 ;
g.self_change(S) = false ;
g.infectable(S) = true ;

% // E
for( i=1:NQUART)
    g.time_dist (E,i) = myQexp( (1.0+i)/NQUART ) / args.sigma;
end
g.state_trans_p(E,1) = args.tau * MAXRAND ; g.state_trans(E,1) = IY ;
g.state_trans_p(E,1) = MAXRAND ;            g.state_trans(E,1) = IA ;

g.state_infect_rate(E) = args.omegaE ;
g.self_change(E) = true ;
g.infectable(E) = false ;

% // IY
rateIY = (1- args.pi) * args.gammaY + args.pi * args.eta;
frac1 =  (1- args.pi) * args.gammaY                        / rateIY ;

for( i=1:NQUART )
    g.time_dist (IY,i) = myQexp( (1.0+i)/NQUART ) / rateIY;
end
g.state_trans_p(IY,1) = frac1 * MAXRAND;  g.state_trans(IY,1) = R ;
g.state_trans_p(IY,1) = MAXRAND;          g.state_trans(IY,1) = IH ;

g.state_infect_rate(IY) = args.omegaY ;
g.self_change(IY) = true ;
g.infectable(IY) = false ;

% // IA
for (i=1:NQUART)
    g.time_dist(IA,i) = myQexp( (1.0+i)/NQUART ) / args.gammaA;
end
g.state_trans_p(IA,1) = MAXRAND;     g.state_trans(IA,1) = R ;

g.state_infect_rate(IA) = args.omegaA ;
g.self_change(IA) = true ;
g.infectable(IA) = false ;

% // IH
rateIH = (1- args.nu) * args.gammaH  +  args.nu * args.mu ;
frac1  = (1- args.nu) * args.gammaH                          / rateIH ;

for( i=1:NQUART)
    g.time_dist(IH,i) = myQexp( (1.0+i)/NQUART ) / rateIH;
end
g.state_trans_p(IH,1) = frac1 * MAXRAND; g.state_trans(IH,1) = R ;
g.state_trans_p(IH,1) = MAXRAND;         g.state_trans(IH,1) = D ;

g.state_infect_rate(IH) = args.omegaH ;
g.self_change(IH) = true ;
g.infectable(IH) = false ;


% // R
g.time_dist(R,:) = NEVER;
g.state_trans_p(R,1) = MAXRAND; g.state_trans(R,1) = R ;

g.state_infect_rate(R) = 0.0 ;
g.self_change(R) = false ;
g.infectable(R) = false ;

% // D
g.time_dist(D,:) = NEVER;
g.state_trans_p(D,1) = MAXRAND ;  g.state_trans(D,1) = D ;

g.state_infect_rate(R) = 0.0 ;
g.self_change(R) = false ;
g.infectable(R) = false ;


end
