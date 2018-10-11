function bLaso = regLasso1(X,Y)
%regLasso1 conduct LASSO method for linear regression, using shooting algorithm (Fu, 1998).
%                                                           (iteration time = 1)
%
%   Positional input:
%     X         A numeric matrix for predictor variables
%     Y         A numeric vector for response variable
%
%   Return values:
%     bLaso     The fitted regression coefficients for each model. 
CountI=1; Errn=1; rrv=1; threv = 0.00002/rrv; Cthrev = 1; p0=0; while  Errn > threv;
W=X; Y0=Y; p = sum( abs(sum(W)) > 0 ); J = sum( abs(sum(W')) > 0 ); 
bpLSE=inv(W'*W)*W'*Y0; sgmL=(Y0-W*bpLSE)'*(Y0-W*bpLSE)/(J-p); kHKB=p*sgmL/(bpLSE'*bpLSE); bpRR=inv(W'*W+kHKB*eye(p))*W'*Y0; beta0=bpRR; if CountI > 1; beta0=bLaso; end; Xr=W; Yr=Y0; Resi=Y0-W*beta0; RSSLSE=Resi'*Resi; Ldaup=100; Vlambda = RSSLSE*[1:100]/200; VDGCV=Vlambda; DBETA=bpRR*ones(1,Ldaup); for ttt=1:Ldaup;
Lda=Vlambda(ttt); RRW=pinv(diag(beta0))/2; IvW=inv(Xr'*Xr+Lda*RRW); PLda=trace( Xr*IvW*Xr' ) ;
betaL=IvW*W'*Yr; ResiL=Y0-W*betaL;  DRS=sum(abs(ResiL)); wgcv=(J*(1-PLda/J)^2);
VDGCV(ttt)=DRS/wgcv; DBETA(:,ttt)=betaL; end; bLaso = DBETA(:,find(VDGCV == min(VDGCV)));  Errn=norm(bLaso(p0+1:p)-beta0(p0+1:p)); CountI=CountI+1; if CountI > Cthrev; Errn=threv; end; end;