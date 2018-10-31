function Y = MyPercentiles( X,p )
%Clone of prctile -- compute percentiles of data in X column by colum
% inputs:
% X -- m x n  data
% p -- np x 1  percentiles
%
% outputs:
% np x n   percentiles for each column of X

m= size(X,1);

p = p(:);

ip = min(max(round(p*m),1),m);

X = sort(X);

Y = X(ip,:);

end

