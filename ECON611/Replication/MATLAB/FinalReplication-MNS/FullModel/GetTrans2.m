function Trans = GetTrans2(q,T)
% Here we use a 3 x 3 transition matrix and increase the 3 -> 1 element and
% the 1 -> 2 element while decreasing 3 -> 2 and 1-> 1 so that we have more
% negative skewness in the transitions while still having the same
% invariant distribution



D = invdistr(T');
T2 = T;
T2(3,1) = T2(3,1) + q;
T2(3,2) = T2(3,2) - q;
T2(1,2) = T2(1,2) + D(3)*q/D(1);
T2(1,1) = T2(1,1) - D(3)*q/D(1);
assert(all(abs(invdistr(T2') - D) < 1e-10))
Trans = T2;

end
