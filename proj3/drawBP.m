function [accepted, t] = drawBP(lambda, t, tau, rho)
d = length(t) - 1;
accepted = false(1,d-1);


%Metropolis-Hasting step
for i = 2:d
    R = rho(i-1)*(t(i+1)-t(i-1));
    t_star =  t(i) - R + 2*R*rand;

    while(t_star < t(i-1) || t_star >  t(i+1))
        t_star =  t(i) - R + 2*R*rand;
    end
    

    num = evalPosteriort(lambda, [t(1:i-1) t_star t(i+1:end)], tau);
    den = evalPosteriort(lambda, t, tau);
    alpha = min(1, num/den);
    U = rand(1);

    if U <= alpha
        t(i) = t_star;
        accepted(i-1) = true;
    end
end



end