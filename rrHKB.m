function bpRR = rrHKB(X,Y)
%rrHKB conduct RR method for linear regression, using Kibria (2003).
%
%   Positional input:
%     X         A numeric matrix for predictor variables
%     Y         A numeric vector for response variable
%
%   Return values:
%     bpRR      The fitted regression coefficients for each model. 
W=X; Y0=Y; p = sum( abs(sum(W)) > 0 ); J = sum( abs(sum(W')) > 0 ); bpLSE=inv(W'*W)*W'*Y0; sgmL=(Y0-W*bpLSE)'*(Y0-W*bpLSE)/(J-p); kHKB=p*sgmL/(bpLSE'*bpLSE); bpRR=inv(W'*W+kHKB*eye(p))*W'*Y0; 