function survBAna2(sDLBCL1,sDLBCL2,sDLBCL3,sDLBCL4)
HRcolor0=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine0=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor0=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine0=[0.18823529779911 0.39215686917305 0.792156875133514]; 
HRcolor=HRcolor0; HRcolorLine=HRcolorLine0; LRcolor=LRcolor0; LRcolorLine=LRcolorLine0; pdd=4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 1: GSE4475 in SPIB-
subplot(7,pdd,[1+pdd*6 1+pdd*4 1+pdd*5]);
Tobs=sDLBCL1(:,1); Delta=sDLBCL1(:,2); Icensor=sDLBCL1(:,3); WCovX=sDLBCL1(:,4:6); nums=1; ppp=0.17; ppS=0.55; 
SCTPs=quantile(WCovX(:,nums),1-ppS);
idxMinus=[ (WCovX(:,nums) <= SCTPs ) ]*1; idxPlus=[ (WCovX(:,nums) > quantile(WCovX(:,nums),ppS) ) ]*1;
SCTP=quantile(WCovX(:,2),ppp);
idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxMinus+idxBAHR == 2); idxLR=find(idxMinus+idxBALR == 2); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(2.1*) BACH2 higher survival better - GSE4475']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); nowDthRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csHRt < XH(jjj)); aaa(bom)=[]; csHRY(aaa)=1-nowDthRt; Achk=[Achk; aaa];
bom=find(csHRt < XH(jjj));  nowDthRt=YH(jjj);
end;
aaa=find(csHRt >= XH(nnn)); Achk=[Achk; aaa]; csHRY(aaa)=1-nowDthRt; AchkH=Achk;
csLRt=sort(LRt( find(LRd == 0) )); csLRY=0*csLRt+1; nnn=sum(XL > -1); nowDtLRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csLRt < XL(jjj)); aaa(bom)=[]; csLRY(aaa)=1-nowDtLRt; Achk=[Achk; aaa];
bom=find(csLRt < XL(jjj));  nowDtLRt=YL(jjj);
end;
aaa=find(csLRt >= XL(nnn)); Achk=[Achk; aaa]; csLRY(aaa)=1-nowDtLRt;
spHRt=0; spHRY=1; nowDthRt=0; nnn=sum(XH > -1); eerr=10^-6; 
for jjj = 1 : nnn;
bfTP=XH(jjj)-eerr; bfY=1-nowDthRt; afTP=XH(jjj); nowDthRt=YH(jjj); afY=1-nowDthRt;
spHRt=[spHRt; bfTP ; afTP]; spHRY=[spHRY; bfY ; afY];
end;
spLRt=0; spLRY=1; nowDtLRt=0; nnn=sum(XL > -1); eerr=10^-6; 
for jjj = 1 : nnn;
bfTP=XL(jjj)-eerr; bfY=1-nowDtLRt; afTP=XL(jjj); nowDtLRt=YL(jjj); afY=1-nowDtLRt;
spLRt=[spLRt; bfTP ; afTP]; spLRY=[spLRY; bfY ; afY];
end;
fsTt=10; LW = 2; PN=plot(0,0,'color', 'w'); hold on; plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;

LV=legend([HRG LRG PN],'BACH2^{-}SPIB^{-} (n = 17)', 'BACH2^{+}SPIB^{-} (n = 53)', 'p = 0.040', 'Location', 'northeast');
LV=legend('boxoff'); set(LV,'FontSize',10); set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]);
ylabel(['Survival  Rate']); xlabel(['Months']); title(['GSE4475'],'FontSize',fsTt, 'FontWeight', 'bold'); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 2: GSE32918 in SPIB-
subplot(7,pdd,[1+pdd*6 1+pdd*4 1+pdd*5]+1);
Tobs=sDLBCL2(:,1); Delta=sDLBCL2(:,2); Icensor=sDLBCL2(:,3); WCovX=sDLBCL2(:,4:6); nums=1; ppp=0.62; ppS=0.75; 
SCTPs=quantile(WCovX(:,nums),1-ppS); idxMinus=[ (WCovX(:,nums) <= SCTPs ) ]*1; idxPlus=[ (WCovX(:,nums) > quantile(WCovX(:,nums),ppS) ) ]*1;
SCTP=quantile(WCovX(:,2),ppp); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxMinus+idxBAHR == 2); idxLR=find(idxMinus+idxBALR == 2); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); 
LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))];  disp(['   ']); 
disp(['(2.2*) BACH2 higher survival better - GSE32918']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); nowDthRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csHRt < XH(jjj)); aaa(bom)=[]; csHRY(aaa)=1-nowDthRt; Achk=[Achk; aaa];
bom=find(csHRt < XH(jjj));  nowDthRt=YH(jjj);
end;
aaa=find(csHRt >= XH(nnn)); Achk=[Achk; aaa]; csHRY(aaa)=1-nowDthRt; AchkH=Achk;
csLRt=sort(LRt( find(LRd == 0) )); csLRY=0*csLRt+1; nnn=sum(XL > -1); nowDtLRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csLRt < XL(jjj)); aaa(bom)=[]; csLRY(aaa)=1-nowDtLRt; Achk=[Achk; aaa];
bom=find(csLRt < XL(jjj));  nowDtLRt=YL(jjj);
end;
aaa=find(csLRt >= XL(nnn)); Achk=[Achk; aaa]; csLRY(aaa)=1-nowDtLRt;
spHRt=0; spHRY=1; nowDthRt=0; nnn=sum(XH > -1); eerr=10^-6; 
for jjj = 1 : nnn;
bfTP=XH(jjj)-eerr; bfY=1-nowDthRt; afTP=XH(jjj); nowDthRt=YH(jjj); afY=1-nowDthRt;
spHRt=[spHRt; bfTP ; afTP]; spHRY=[spHRY; bfY ; afY];
end;
spLRt=0; spLRY=1; nowDtLRt=0; nnn=sum(XL > -1); eerr=10^-6; 
for jjj = 1 : nnn;
bfTP=XL(jjj)-eerr; bfY=1-nowDtLRt; afTP=XL(jjj); nowDtLRt=YL(jjj); afY=1-nowDtLRt;
spLRt=[spLRt; bfTP ; afTP]; spLRY=[spLRY; bfY ; afY];
end;
fsTt=10; LW = 2; PN=plot(0,0,'color', 'w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;

LV=legend([HRG LRG PN],'BACH2^{-}SPIB^{-} (n = 28)', 'BACH2^{+}SPIB^{-} (n = 14)', 'p = 0.155', 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',10); set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]);
ylabel(['Survival  Rate']); xlabel(['Months']); title(['GSE32918'],'FontSize',fsTt, 'FontWeight', 'bold');  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 3: GSE31312 in SPIB-
subplot(7,pdd,[1 1+pdd 1+pdd*2]);
Tobs=sDLBCL3(:,1); Delta=sDLBCL3(:,2); Icensor=sDLBCL3(:,3); WCovX=sDLBCL3(:,4:6); nums=1; ppp=0.21; ppS=0.5; 
SCTPs=quantile(WCovX(:,nums),1-ppS); idxMinus=[ (WCovX(:,nums) <= SCTPs ) ]*1; idxPlus=[ (WCovX(:,nums) > quantile(WCovX(:,nums),ppS) ) ]*1;
SCTP=quantile(WCovX(:,2),ppp); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxMinus+idxBAHR == 2); idxLR=find(idxMinus+idxBALR == 2); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(2.3*) BACH2 higher survival better - GSE31312']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']); 

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); nowDthRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csHRt < XH(jjj)); aaa(bom)=[]; csHRY(aaa)=1-nowDthRt; Achk=[Achk; aaa];
bom=find(csHRt < XH(jjj));  nowDthRt=YH(jjj);
end;
aaa=find(csHRt >= XH(nnn)); Achk=[Achk; aaa]; csHRY(aaa)=1-nowDthRt; AchkH=Achk;
csLRt=sort(LRt( find(LRd == 0) )); csLRY=0*csLRt+1; nnn=sum(XL > -1); 
nowDtLRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csLRt < XL(jjj)); aaa(bom)=[]; csLRY(aaa)=1-nowDtLRt; Achk=[Achk; aaa];
bom=find(csLRt < XL(jjj));  nowDtLRt=YL(jjj);
end;
aaa=find(csLRt >= XL(nnn)); Achk=[Achk; aaa]; csLRY(aaa)=1-nowDtLRt; spHRt=0; spHRY=1; nowDthRt=0; nnn=sum(XH > -1); eerr=10^-6; 
for jjj = 1 : nnn;
bfTP=XH(jjj)-eerr; bfY=1-nowDthRt; afTP=XH(jjj); nowDthRt=YH(jjj); afY=1-nowDthRt;
spHRt=[spHRt; bfTP ; afTP]; spHRY=[spHRY; bfY ; afY];
end;
spLRt=0; spLRY=1; nowDtLRt=0; nnn=sum(XL > -1); eerr=10^-6;
for jjj = 1 : nnn;
bfTP=XL(jjj)-eerr; bfY=1-nowDtLRt; afTP=XL(jjj); nowDtLRt=YL(jjj); afY=1-nowDtLRt;
spLRt=[spLRt; bfTP ; afTP]; spLRY=[spLRY; bfY ; afY];
end;
fsTt=10; LW = 2; PN=plot(0,0,'color', 'w'); hold on; plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;

LV=legend([HRG LRG PN],'BACH2^{-}SPIB^{-} (n = 61)', 'BACH2^{+}SPIB^{-} (n = 174)', 'p = 0.0275', 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',10); set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]);
ylabel(['Survival  Rate']); xlabel(['Months']); title(['GSE31312'],'FontSize',fsTt, 'FontWeight', 'bold');  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 4: GSE10846 in SPIB-
subplot(7,pdd,[ [1 1+pdd 1+pdd*2]+1 ]);
Tobs=sDLBCL4(:,1); Delta=sDLBCL4(:,2); Icensor=sDLBCL4(:,3); WCovX=sDLBCL4(:,4:6); nums=1; ppp=0.21; ppS=0.54; 
SCTPs=quantile(WCovX(:,nums),1-ppS); idxMinus=[ (WCovX(:,nums) <= SCTPs ) ]*1; idxPlus=[ (WCovX(:,nums) > quantile(WCovX(:,nums),ppS) ) ]*1;
SCTP=quantile(WCovX(:,2),ppp); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxMinus+idxBAHR == 2); idxLR=find(idxMinus+idxBALR == 2); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); 
LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(2.4*) BACH2 higher survival better - GSE10846']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']); 

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); nowDthRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csHRt < XH(jjj)); aaa(bom)=[]; csHRY(aaa)=1-nowDthRt; Achk=[Achk; aaa];
bom=find(csHRt < XH(jjj));  nowDthRt=YH(jjj);
end;
aaa=find(csHRt >= XH(nnn)); Achk=[Achk; aaa]; csHRY(aaa)=1-nowDthRt; AchkH=Achk;
csLRt=sort(LRt( find(LRd == 0) )); csLRY=0*csLRt+1; nnn=sum(XL > -1); nowDtLRt=0; bom=[]; Achk=[];
for jjj = 1 : nnn;
aaa=find(csLRt < XL(jjj)); aaa(bom)=[]; csLRY(aaa)=1-nowDtLRt; Achk=[Achk; aaa];
bom=find(csLRt < XL(jjj));  nowDtLRt=YL(jjj);
end;
aaa=find(csLRt >= XL(nnn)); Achk=[Achk; aaa]; csLRY(aaa)=1-nowDtLRt;
spHRt=0; spHRY=1; nowDthRt=0; nnn=sum(XH > -1); eerr=10^-6;
for jjj = 1 : nnn;
bfTP=XH(jjj)-eerr; bfY=1-nowDthRt; afTP=XH(jjj); nowDthRt=YH(jjj); afY=1-nowDthRt;
spHRt=[spHRt; bfTP ; afTP]; spHRY=[spHRY; bfY ; afY];
end;
spLRt=0; spLRY=1; nowDtLRt=0; nnn=sum(XL > -1); eerr=10^-6;
for jjj = 1 : nnn;
bfTP=XL(jjj)-eerr; bfY=1-nowDtLRt; afTP=XL(jjj); nowDtLRt=YL(jjj); afY=1-nowDtLRt;
spLRt=[spLRt; bfTP ; afTP]; spLRY=[spLRY; bfY ; afY];
end;
fsTt=10; LW = 2; PN=plot(0,0,'color', 'w'); hold on; plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;

LV=legend([HRG LRG PN],'BACH2^{-}SPIB^{-} (n = 39)', 'BACH2^{+}SPIB^{-} (n = 151)', 'p = 0.0025', 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',10); set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]);
ylabel(['Survival  Rate']); xlabel(['Months']); title(['GSE10846'],'FontSize',fsTt, 'FontWeight', 'bold');  


