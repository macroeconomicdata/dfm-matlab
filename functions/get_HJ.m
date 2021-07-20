%--------------------------------------------------------------------------

% Helper matrix J for observation equation

function HJ = get_HJ(H,frq,is_diff,p)
[k,m] = size(H);
lags = frq;
lags(is_diff,:) = arrayfun(@(x)(2*x-1),frq(is_diff,:));
pp = max([lags;p]);
HJ = zeros(k, m*pp);
HJ(:,1:m) = H; %for high frequency data
idx = find(frq>1); %for low frequency data
for i = 1:size(idx,1)
    j = idx(i);
    if is_diff(j)
        HJ(j,1:m*(2*frq(j)-1)) = H(j,:)*kron([1:frq(j),(frq(j)-1):-1:1]/frq(j), eye(m));
    else
        HJ(j,1:m*frq(j)) = H(j,:)*kron(ones(1,frq(j))/frq(j),eye(m));
    end
end 
return
end