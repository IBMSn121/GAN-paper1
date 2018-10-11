function survTAna3(sOSBC3a,sOSBC3b)
pdd=4; fsTtF=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MM1 - GSE9782
subplot(7,pdd,[1 1+pdd 1+pdd*2]);
Tobs=sOSBC3a(:,1); Delta=sOSBC3a(:,2); Icensor=sOSBC3a(:,3); WCovX=sOSBC3a(:,4:6);
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg;
qcp=0.83; SCTP=quantile(ScorePH,qcp); 
HRd=Delta(find( (ScorePH > SCTP) ) ); HRt=Tobs(find( (ScorePH > SCTP) ) ); 
LRd=Delta(find( (ScorePH <= SCTP) ) ); LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(6.1*) Trio''s prognosis for MM - GSE9782']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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

LV=legend([HRG LRG PN],{'high risk (n = 45)', 'low risk (n = 219)', 'p = 3.32E-06'}, 'Location', 'southwest'); LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE9782 - MM'],'FontSize',fsTt, 'FontWeight', 'bold');   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MM2 - GSE24080
subplot(7,pdd,[ [1 1+pdd 1+pdd*2]+1 ]);
Tobs=sOSBC3b(:,1); Delta=sOSBC3b(:,2); Icensor=sOSBC3b(:,3); WCovX=sOSBC3b(:,4:6);
bphreg = coxphfit(WCovX, Tobs, 'Censoring', Icensor); ScorePH=WCovX*bphreg; 
qcp=0.3; SCTP=quantile(ScorePH,qcp); HRd=Delta(find( (ScorePH > SCTP) ) ); 
HRt=Tobs(find( (ScorePH > SCTP) ) ); LRd=Delta(find( (ScorePH <= SCTP) ) );
LRt=Tobs(find( (ScorePH <= SCTP) ) ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1); ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(6.2*) Trio''s prognosis for MM - GSE24080']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

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

LV=legend([HRG LRG PN],{'high risk (n = 391)', 'low risk (n = 168)', 'p = 0.011'}, 'Location', 'southwest'); LV=legend('boxoff'); set(LV,'FontSize',11); set(gcf,'Color',[1,1,1]); 
axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE24080 - MM'],'FontSize',fsTt, 'FontWeight', 'bold'); 




