function output_value = evalPosteriort(lambda, t, tau)
d = length(t)-1;


n_save = zeros(d, 1);
for i = 1:d
    num_disasters = find(tau >= t(i) & tau<t(i+1));
    n_save(i) = length(num_disasters);
end

t_diff = zeros(d, 1);
for l = 1:d
    t_diff(l) = t(l+1) - t(l);
end


output_value = exp(sum(log(lambda).*n_save + log(t_diff) - lambda.*t_diff));
end