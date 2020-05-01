

%%
function g = seir_init(nS, NSTATES, g)
% 	unsigned int i, source;
S=0;
% // g.t = 0.0;

for(i=1:NSTATES)
    g.s(i) = 0 ;
end

% // initialize
for (i = 1: g.n)
    nd(i).heap = NONE;
    nd(i).time = DBL_MAX; % // to a large value
    nd(i).state = nd(i).state % g.ngroup ; % // reset to modulo ngroup
    g.s( nd(i).state ) = g.s( nd(i).state )+1;
end

for (i=1:nS)
    % // get & infect the source
    source = pcg_32_bounded(g.n);
    
    while( nd(source).time==0.0 )
        source = pcg_32_bounded(g.n);
    end;
    nd(source).time = 0.0;
    nd(source).next_state = nd(source).state+ 1* g.ngroup ; % // assuming for each group, 0 is S, 1 is E.
    nd(source).heap = i;
    g.heap(g.nheap == i ) = source;
end

end
