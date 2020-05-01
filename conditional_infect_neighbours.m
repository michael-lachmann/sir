


function conditional_infect_neighbours(t_me, me,  rate)
% // go through my neighbors
for (i = 1:nd(me).deg)
    you = nd(me).nb(i);
    if (g.infectable(nd(you).state ))  % // if you can be infected
        t = g.now + ...
            g.time_dist( nd(you).state ,RANDOM_QUART()) / (g.weight( nd(me).w(i)  ) *rate ) ; % // get the infection time
        if ((t < t_me ) &&   ...  % // When me becomes resistant
                (t < nd(you).time))     % // when you become exposed, if you do
            % 				DEBUG( printf("%g %g inf:%d\n",t_me,t,you) ) ;
            nd(you).time = t;
            nd(you).next_state = next_state( nd(you).state ) ; % // ***** NOT RANDOM YET, because we always switch to E.
            if (nd(you).heap == NONE)  % // if not listed before, then extend the heap
                g.nheap = g.nheap+1;
                g.heap(g.nheap) = you;
                nd(you).heap = g.nheap;
                
            end
            up_heap(nd(you).heap); % // this works bcoz the only heap relationship
            % // that can be violated is the one between you and its parent
        end
    end
end
end

