function survTAna1(sDLBCL1,sDLBCL2,sDLBCL3,sDLBCL4)
pdd=4; fsTtF=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 1: GSE4475
subplot(7,pdd,[1 1+pdd 1+pdd*2]);
Tobs=sDLBCL1(:,1); Delta=sDLBCL1(:,2); Icensor=sDLBCL1(:,3); WCovX=sDLBCL1(:,4:6);
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg;
qcp=0.73; SCTP=quantile(ScorePH,qcp); 
HRd=Delta(find( (ScorePH > SCTP) ) ); HRt=Tobs(find( (ScorePH > SCTP) ) ); 
LRd=Delta(find( (ScorePH <= SCTP) ) ); LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(4.1*) Trio''s prognosis for DLBCL - GSE4475']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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
HRcolor=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine=[0.18823529779911 0.39215686917305 0.792156875133514];   
fsTt=fsTtF; LW = 2; PN=plot(0,0,'color', 'w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine);
LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine);
hold off;

LV=legend([HRG LRG PN],{'high risk (n = 42)', 'low risk (n = 113)', 'p = 5.44E-06'}, 'Location', 'southwest'); LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE4475'],'FontSize',fsTt);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 2: GSE32918 in SPIB
subplot(7,pdd,[ [1 1+pdd 1+pdd*2]+1 ]);
Tobs=sDLBCL2(:,1); Delta=sDLBCL2(:,2); Icensor=sDLBCL2(:,3); WCovX=sDLBCL2(:,4:6);
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg; 
qcp=0.15; SCTP=quantile(ScorePH,qcp); HRd=Delta(find( (ScorePH > SCTP) ) ); 
HRt=Tobs(find( (ScorePH > SCTP) ) ); LRd=Delta(find( (ScorePH <= SCTP) ) );
LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(4.2*) Trio''s prognosis for DLBCL - GSE32918']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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
HRcolor=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine=[0.18823529779911 0.39215686917305 0.792156875133514];   
fsTt=fsTtF; LW = 2; PN=plot(0,0,'color', 'w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine);
LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine);
hold off;

LV=legend([HRG LRG PN],{'high risk (n = 142)', 'low risk (n = 25)', 'p = 0.033'}, 'Location','southwest');
LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE32918'],'FontSize',fsTt); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 3: GSE31312 in SPIB
subplot(7,pdd,[1+pdd*6 1+pdd*4 1+pdd*5]);
Tobs=sDLBCL3(:,1); Delta=sDLBCL3(:,2); Icensor=sDLBCL3(:,3); WCovX=sDLBCL3(:,4:6); 
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg;
qcp=0.5; SCTP=quantile(ScorePH,qcp); HRd=Delta(find( (ScorePH > SCTP) ) );
HRt=Tobs(find( (ScorePH > SCTP) ) ); LRd=Delta(find( (ScorePH <= SCTP) ) ); 
LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); 
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(4.3*) Trio''s prognosis for DLBCL - GSE31312']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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
HRcolor=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine=[0.18823529779911 0.39215686917305 0.792156875133514];   
fsTt=fsTtF; LW = 2; PN=plot(0,0,'color', 'w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine);
LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine);
hold off;

LV=legend([HRG LRG PN],{'high risk (n = 235)', 'low risk (n = 235)', 'p = 0.0007'}, 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE31312'],'FontSize',fsTt);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%panel 4: GSE10846 in SPIB
subplot(7,pdd,[1+pdd*6 1+pdd*4 1+pdd*5]+1);
Tobs=sDLBCL4(:,1); Delta=sDLBCL4(:,2); Icensor=sDLBCL4(:,3); WCovX=sDLBCL4(:,4:6);
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg;
qcp=0.78; SCTP=quantile(ScorePH,qcp); 
HRd=Delta(find( (ScorePH > SCTP) ) ); HRt=Tobs(find( (ScorePH > SCTP) ) ); 
LRd=Delta(find( (ScorePH <= SCTP) ) ); LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(4.4*) Trio''s prognosis for DLBCL - GSE10846']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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
HRcolor=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine=[0.18823529779911 0.39215686917305 0.792156875133514];   
fsTt=fsTtF; LW = 2; PN=plot(0,0,'color', 'w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine);
LRG=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine);
hold off;

LV=legend([HRG LRG PN],{'high risk (n = 91)', 'low risk (n = 321)', 'p = 0.0003'}, 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE10846'],'FontSize',fsTt);  



