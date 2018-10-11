function [Beta,Pv] = regAWTE(X,Y)
%regAWTE conduct AWTE method for linear regression.
%
%   Positional input:
%     X         A numeric matrix for predictor variables (the last column: common interaction component)
%     Y         A numeric vector for response variable
%
%   Return values:
%     Beta             The fitted regression coefficients for each model. 
%     Pv               The p-values corresponding to Beta using T-tests.

W=X; Y0=Y; p = sum( abs(sum(W)) > 0 ); J = sum( abs(sum(W')) > 0 ); 
AWTE=zeros(p,1); a0h=zeros(p,1); terror=0.00000001; G=zeros(p,J); I_p=eye(p);
X0=W(:,p); PupI=find(X0 >= median(X0)+terror); PlwI=find(X0 < median(X0)+terror); Yup=mean(Y0(PupI)); Ylw=mean(Y0(PlwI)); Xup=mean(X0(PupI)); Xlw=mean(X0(PlwI)); cwp=(Xup-Xlw); Kstar=(Yup-Ylw)/cwp; 
gup=(X0 >= median(X0)+terror)/sum(X0 >= median(X0)+terror); glw=(X0 < median(X0)+terror)/sum(X0 < median(X0)+terror); G(p,:)=gup-glw; for t = 1 : p-1; H0=W(:,t); a0=[mean(H0(PupI))-mean(H0(PlwI))]/cwp; X0=H0-a0*W(:,p); a0h(t)=a0; upI=find(X0 >= median(X0)+terror); lwI=find(X0 < median(X0)+terror); Yup=mean(Y0(upI)); Ylw=mean(Y0(lwI)); Xup=mean(X0(upI)); Xlw=mean(X0(lwI)); AWTE(t)=(Yup-Ylw)/(Xup-Xlw); gup=(X0 >= median(X0)+terror)/sum(X0 >= median(X0)+terror); glw=(X0 < median(X0)+terror)/sum(X0 < median(X0)+terror); G(t,:)=gup-glw; end; AWTE(p)=Kstar-a0h'*AWTE; Wa=G*(W-W(:,p)*a0h'); DW=diag(diag(Wa)); Ia=I_p-[zeros(p-1,p); a0h']; AW_m=Ia*inv(DW)*G; Is_AW=eye(J)-AW_m'*inv(AW_m*AW_m')*AW_m; GMSE_AW=Y0'*Is_AW*Y0/(J-p); Sgm2AW=diag(  AW_m*AW_m'   )'; NT_AWTE=(AWTE')./sqrt(Sgm2AW*GMSE_AW); NP_AWTE=2*(1-tcdf(abs(NT_AWTE),J-p)); Beta=AWTE; Pv=NP_AWTE';
