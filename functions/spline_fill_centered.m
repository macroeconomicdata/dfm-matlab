function y = spline_fill_centered(x)
mx = mean(x, 'omitnan');
T = length(x);
ind = 1:T;
ind_obs = ind(~isnan(x)); % index of all observations
fst = ind_obs(1); % index of first observation
lst = ind_obs(length(ind_obs)); % index of last observation
y = mx*ones(T,1); % if at beginning or end just use the mean
y(fst:lst) = spline(ind_obs,x(ind_obs),fst:lst); % spline from first to 
% last observation, all others are mean
end