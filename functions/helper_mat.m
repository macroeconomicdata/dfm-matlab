%--------------------------------------------------------------------------

% Helper matrix J for observation equation

function J = helper_mat(fq,isdiff,r,m)
J = zeros(r,m);
    if isdiff
        J(:,1:r*(2*fq-1)) = kron([1:fq,(fq-1):-1:1]/fq, eye(r));
    else
        J(:,1:r*fq) = kron(ones(1,fq)/fq,eye(r));
    end 
return
end