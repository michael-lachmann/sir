

function res = next_state(s, g, MAXRAND)
i=0;
while g.state_trans_p(s,i) < MAXRAND
    if rand < g.state_trans_p(s,i)
        res = g.state_trans(s,i);
        break
    end
    i = i+1 ;
end
res =  g.state_trans(s,i) ;
end