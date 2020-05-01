


function conditional_infect_neighbours = reinfect_s_I_neighbours(g, conditional_infect_neighbours)
% 	unsigned int i  ;
for i = 1:g.n
    if ( g.state_infect_rate( nd(i).state ) > 0)
        conditional_infect_neighbours(  nd(i).time, i,  g.state_infect_rate( nd(i).state )) ;
    end
end
end