function y = spline_fill_plot(x)
T = length(x);
ind = 1:T;
ind_obs = ind(~isnan(x));
fst = ind_obs(1);
lst = ind_obs(length(ind_obs));
y = NaN(T,1);
y(fst:lst) = spline(ind_obs,x(ind_obs),fst:lst);
end