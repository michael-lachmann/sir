

function epi_timestep()
% unsigned int me ;
% state next_s, s ;
% next_t ;

me = g.heap(1) ;
g.now = nd(me).time ;

% // Advance node to its next state
% //		nd(me).time  = nd(me).next_time ;

g.s( nd(me).state      ) = g.s( nd(me).state      )-1;
g.s( nd(me).next_state ) = g.s( nd(me).next_state )+1;
% 		DEBUG( printf("%d ",me) ) ;
% 		DEBUG( print_state( nd(me).state) ) ;
% 		DEBUG( printf("->") ) ;
% 		DEBUG( print_state( nd(me).next_state) ) ;
% 		DEBUG( printf("\n") ) ;
nd(me).state = nd(me).next_state;


% // state detemines what happens next
s = nd(me).state ;

% // Node just switched to state s.
% // What to do next with it?
% // If it is in a state that will change in the future on its own, then schedule that.
% // If it is a state that can infect others, schedule those infections.
if( g.self_change(s) ) % // State is one that will change on its own in the future
    % // time of next own change
    %     next_t = g.now + g.time_dist(s)[RANDOM_QUART()] ;
    next_t = g.now + g.time_dist(s, randi( size(g.time_dist,2)) ) ;
    % // choose next state
    next_s = next_state(s) ;
    
    % // If infection rate is bigger than 0, test for infections.
    if( g.state_infect_rate(s) > EPSILON )
        conditional_infect_neighbours( next_t, me, g.state_infect_rate(s)) ; % // only infect neighbours that get infected before switch to next state.
    end
    % // schedule next event for node
    nd(me).next_state = next_s ;
    nd(me).time       = next_t ;
    down_heap( nd(me).heap ) ;
else  % // State doesn't change on its own, schedule a change for NEVER.
    % // We switched to a state that should never change
    nd(me).time = NEVER ;
    % // next_state stays what it is
    down_heap( nd(me).heap ) ;
end
end
