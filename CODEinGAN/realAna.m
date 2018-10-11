function realAna(dX,dY)
TrueBBPPVC4_L2L6L=[1 1 1 1 1 1 1 0 1 1; 0 0 0 0 0 1 1 0 0 0; 0 0 0 1 0 0 1 0 0 0; 0 1 0 1 1 1 1 1 0 0; 0 0 0 0 0 0 0 1 0 0]'; idxPPVC4=[1 2 3 5 4]; TrueBBO1_5NA=TrueBBPPVC4_L2L6L(:,idxPPVC4); idx1=[2 3 4 5 1]; ns=133; nga=5; ngb=10; p=nga; nnn=ngb; nel=nnn; J=ns; indXa=idx1; TrueBB=TrueBBO1_5NA(:,idx1); LPCYo=dY; W=dX; p0=0; warning off;

%Inference
Bestmt_AWTE=ones(nel,p); Bpvalu_AWTE=ones(nel,p); Bsigsc_AWTE=ones(nel,p);
Bestmt_LASSO=ones(nel,p); Bpvalu_LASSO=ones(nel,p); Bsigsc_LASSO=ones(nel,p); 
Bestmt_RRE=ones(nel,p); Bpvalu_RRE=ones(nel,p); Bsigsc_RRE=ones(nel,p); 
Bestmt_LSE=ones(nel,p); Bpvalu_LSE=ones(nel,p); Bsigsc_LSE=ones(nel,p);  
for ggg=1:nnn; Y0=LPCYo(:,ggg); bpLSE=inv(W'*W)*W'*Y0; BpLaso=regLasso1(W, Y0); bpRR=rrHKB(W, Y0);
[Beta,Pv] = regAWTE(W, Y0); AWTE=Beta; sgmL=(Y0-W*bpLSE)'*(Y0-W*bpLSE)/(J-p); kHKB=p*sgmL/(bpLSE'*bpLSE); RR_m=inv(W'*W+kHKB*eye(p+p0))*W';Is_RR=eye(J)-RR_m'*inv(RR_m*RR_m')*RR_m; GMSE_RR=Y0'*Is_RR*Y0/(J-p); Sgm2RR=diag(  RR_m*RR_m'   )'; MSE_LS=(Y0-W*bpLSE)'*(Y0-W*bpLSE)/(J-p); Sigma2=diag(inv(W'*W))'; MSE_LA=(Y0-W*BpLaso)'*(Y0-W*BpLaso)/(J-p); T_LASSO=(BpLaso')./sqrt(Sigma2*MSE_LA); NT_RRE=(bpRR')./sqrt(Sgm2RR*GMSE_RR); T_LSE=(bpLSE')./sqrt(Sigma2*MSE_LS); NP_AWTE=Pv; P_LASSO=2*(1-tcdf(abs(T_LASSO),J-p)); NP_RRE=2*(1-tcdf(abs(NT_RRE),J-p)); P_LSE=2*(1-tcdf(abs(T_LSE),J-p)); Bestmt_AWTE(ggg,:)=AWTE'; Bpvalu_AWTE(ggg,:)=NP_AWTE; Bsigsc_AWTE(ggg,:)=(NP_AWTE < 0.05); Bestmt_LASSO(ggg,:)=BpLaso'; Bpvalu_LASSO(ggg,:)=P_LASSO; Bsigsc_LASSO(ggg,:)=(P_LASSO < 0.05);  Bestmt_RRE(ggg,:)=bpRR'; Bpvalu_RRE(ggg,:)=NP_RRE; Bsigsc_RRE(ggg,:)=(NP_RRE < 0.05); Bestmt_LSE(ggg,:)=bpLSE'; Bpvalu_LSE(ggg,:)=P_LSE; Bsigsc_LSE(ggg,:)=(P_LSE < 0.05);  
end; 
%AUC calculations
TAnum=50; Iref=reshape(TrueBB,TAnum,1); 
IAW=reshape(abs(Bestmt_AWTE),TAnum,1); [X1,Y1,T1,AUC1,OPTROCPT1]=perfcurve(Iref,IAW,1); [AUC1 OPTROCPT1];
ILS=reshape(abs(Bestmt_LSE),TAnum,1); [X4,Y4,T4,AUC4,OPTROCPT4]=perfcurve(Iref,ILS,1); [AUC4 OPTROCPT4];
IRR=reshape(abs(Bestmt_RRE),TAnum,1); [X3,Y3,T3,AUC3,OPTROCPT3]=perfcurve(Iref,IRR,1); [AUC3 OPTROCPT3];
ILA=reshape(abs(Bestmt_LASSO),TAnum,1); [X2,Y2,T2,AUC2,OPTROCPT2]=perfcurve(Iref,ILA,1); [AUC2 OPTROCPT2];
AUC=[AUC1 AUC2 AUC3 AUC4]'; 
disp(' ')
    disp('(1*) Summary of the AUCs of the four different methods')
	fprintf('-----------------------------------------\n');
    disp('AWTE        LASSO        RRE        LSE ' );
    fprintf('-----------------------------------------\n');
    fprintf('%10s %10.0f %10.0f %10.0f %10.2f\n', num2str( round(AUC'+0.0001,2) ) ); 
    fprintf(' \n');
	fprintf('-----------------------------------------\n');

%plot
ck=48; cc=[0:1/ck:1]'; cc=[cc cc cc]; cclr2=cc(41,:); cclr1=cclr2; cc=hsv(ck); cclr0=cc(44,:);  cclr1(2)=0.65; cclr3=cclr2; cclr3(1)=0.75; ssk=5; aad=2; fsTt=11; fsTtA=8; uby=1.1; lby=-0.1; aa1=1; aa2=3; bb1=4; bb2=6; cclr=[0 0 0]; C1=([1:11]-1)/10; C2=C1; h=area(C1,C2); rrcc=0.5; h.FaceColor = [1 1 1]*rrcc; h.EdgeColor = [1 1 1]*rrcc; hold on;
AW=plot(X1,Y1,'-','color', cclr0,'LineWidth', 1.2,'MarkerSize', 6);
LA=plot(X2,Y2,'--','color', cclr1, 'LineWidth', 1,'MarkerSize', 6);
RR=plot(X3,Y3,':','color', cclr2, 'LineWidth', 1.2,'MarkerSize', 6);
LS=plot(X4,Y4,'-.','color', cclr3, 'LineWidth', 1,'MarkerSize', 6);
hold off;
legend([AW LA RR LS], {'AWTE', 'LASSO', 'RRE', 'LSE'}, 'Location','southeast','FontSize',9);
legend('boxoff'); set(gcf,'Color',[1,1,1]); axis tight;
xlabel('1 - Specificity','FontSize',fsTt, 'FontWeight', 'bold');
ylabel('Sensitivity','FontSize',fsTt, 'FontWeight', 'bold');


%Bestmt_AWTE'; Bpvalu_AWTE';
ABeta=Bestmt_AWTE(:,[5 1:4])'; APvle=Bpvalu_AWTE(:,[5 1:4])';
%BCL6	SPIB	BACH2	IRF8	OCT2
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
disp(' ')
    disp('(2*) AWTE-induced GAN inference - association parameter estimations')
	fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
disp(['        p21         MYC         p53        BCL2        NFKB1       IRF4       Blimp1        AID        p27         ATR']);
    fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
disp(['BCL6   ' num2str( round(ABeta(1,:)+0.0001,3) )])
disp(['SPIB   ' num2str( round(ABeta(2,:)+0.0001,3) )])
disp(['BACH2 ' num2str( round(ABeta(3,:)+0.0001,3) )])
disp(['IRF8  ' num2str( round(ABeta(4,:)+0.0001,3) )])
disp(['OCT2  ' num2str( round(ABeta(5,:)+0.0001,3) )])
fprintf('----------------------------------------------------------------------------------------------------------------------------\n');


disp(' ')
    disp('(3*) AWTE-induced GAN inference - Statistical testing P-values')
	fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
disp(['        p21         MYC         p53        BCL2        NFKB1       IRF4       Blimp1        AID        p27         ATR']);
    fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
disp(['BCL6   ' num2str( round(APvle(1,:)+0.0001,3) )])
disp(['SPIB   ' num2str( round(APvle(2,:)+0.0001,3) )])
disp(['BACH2  ' num2str( round(APvle(3,:)+0.0001,3) )])
disp(['IRF8   ' num2str( round(APvle(4,:)+0.0001,3) )])
disp(['OCT2   ' num2str( round(APvle(5,:)+0.0001,3) )])
fprintf('----------------------------------------------------------------------------------------------------------------------------\n');