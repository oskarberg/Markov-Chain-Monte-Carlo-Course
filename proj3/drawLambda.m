function lambda = posteriorLambda(theta, t, tau)
    tDiff = t(2:end) - t(1:end-1);
    d = length(t) - 1;

    n = zeros(1, d);
    for i = 1:d
        n(i) = sum((t(i) <= tau) & (tau < t(i+1)));
    end
    lambda = gamrnd(n' + 2, 1./(theta + tDiff'));

end
