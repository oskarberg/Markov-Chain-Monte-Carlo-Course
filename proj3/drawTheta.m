function theta = drawTheta(lambda, psi)

    d = length(lambda);
    theta = gamrnd(2*d + 2, 1./(psi + sum(lambda)));
end