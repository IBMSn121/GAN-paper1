function survBAna1(sDLBCL1,sDLBCL2,sDLBCL3,sDLBCL4)
HRcolor0=[0.874509811401367 0.10196078568697 0.968627452850342]; HRcolorLine0=[0.874509811401367 0.133333340287209 0.34901961684227];  
LRcolor0=[0.254901975393295 0.643137276172638 0.901960790157318]; LRcolorLine0=[0.18823529779911 0.39215686917305 0.792156875133514]; 
HRcolorLine1=[0.87058824300766 0.658823549747467 0.14509804546833]; HRcolor1=[0.796078443527222 0.874509811401367 0.164705887436867];  
LRcolorLine1=[0.141176477074623 0.623529434204102 0.384313732385635]; LRcolor1=[0.23137255012989 0.831372559070587 0.109803922474384]; 
HRcolor1a=[0.980392158031464 0.839215695858002 0.0666666701436043];
HRcolorLine2=[0.627451002597809 0.0901960805058479 0.764705896377563];  HRcolor2=[0.996078431606293 0.450980395078659 0.996078431606293];
LRcolorLine2=[0 0.329411774873734 0.815686285495758]; LRcolor2=[0.0392156876623631 0.843137264251709 0.843137264251709]; 
LRcolorLine2a=[0 0.447058826684952 0.74117648601532];
HRcolor1=HRcolor0; HRcolorLine1=HRcolorLine0;  LRcolor1=LRcolor0; LRcolorLine1=LRcolorLine0; 
HRcolor2=[0.627451002597809 0.0901960805058479 0.764705896377563]; HRcolorLine2=[0.996078431606293 0.450980395078659 0.996078431606293]; LRcolor2=[0 0.329411774873734 0.815686285495758]; LRcolorLine2=[0.0392156876623631 0.843137264251709 0.843137264251709]; pdd=6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%left-panel: GSE4475-GSE31312
subplot(4,pdd,[1 1+pdd 1+pdd*2 [1 1+pdd 1+pdd*2]+1]);
%%%whole GSE4475; for BACH2
%%%
Tobs=sDLBCL1(:,1); Delta=sDLBCL1(:,2); Icensor=sDLBCL1(:,3); WCovX=sDLBCL1(:,4:6);
SCTP=quantile(WCovX(:,2),0.8); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxBAHR == 1); idxLR=find(idxBALR == 1); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(1.1*) Controversy of BACH2''s role - GSE4475']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); 
nowDthRt=0; bom=[]; Achk=[];
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
HRcolor=HRcolor2; HRcolorLine=HRcolorLine2;  
LRcolor=LRcolor2; LRcolorLine=LRcolorLine2;
fsTt=10; LW = 2; PN1=plot(-0.5,-0.5,'color','w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG1=plot(spHRt,spHRY,':','LineWidth',LW,'Color',HRcolorLine);
LRG1=plot(spLRt,spLRY,':','LineWidth',LW,'Color',LRcolorLine);

%%%whole GSE31312; for BACH2
%%%
Tobs=sDLBCL3(:,1); Delta=sDLBCL3(:,2); Icensor=sDLBCL3(:,3); WCovX=sDLBCL3(:,4:6);
SCTP=quantile(WCovX(:,2),0.15); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxBAHR == 1); idxLR=find(idxBALR == 1); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2));
YH=fH; XH=xH; XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']); 
disp(['(1.2*) Controversy of BACH2''s role - GSE31312']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); 
nowDthRt=0; bom=[]; Achk=[];
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
HRcolor=HRcolor1; HRcolorLine=HRcolorLine1;  
LRcolor=LRcolor1; LRcolorLine=LRcolorLine1;  
fsTt=10; LW = 2; 
PN2=plot(-0.6,-0.6, '.', 'color', 'w'); PNc=plot(-0.55,-0.55, '.', 'color', 'w');
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG2=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG2=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;
LV=legend([HRG2 LRG2 PNc PN2 HRG1 LRG1], {'BACH2^{-}  (n = 71)', 'BACH2^{+} (n = 399)',  '      ', '      ', 'BACH2^{-} (n = 124)', 'BACH2^{+} (n = 31)'}, 'Location', 'southwest');
LV=legend('boxoff'); set(LV,'FontSize',11);
text(2.5,0.35, 'GSE31312 (p = 0.0385)','FontWeight','bold','FontSize',11,'HorizontalAlignment','left'); 
text(2.5,0.16, 'GSE4475 (p = 0.303)','FontWeight','bold','FontSize',11,'HorizontalAlignment','left'); 
set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]); ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE4475 - GSE31312'],'FontSize',fsTt, 'FontWeight', 'bold');  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%right-panel: GSE32918-GSE10846
subplot(4,pdd,[ [1 1+pdd 1+pdd*2]+2 [1 1+pdd 1+pdd*2]+3 ]);
%%%whole GSE32918; for BACH2
%%%
Tobs=sDLBCL2(:,1); Delta=sDLBCL2(:,2); Icensor=sDLBCL2(:,3); WCovX=sDLBCL2(:,4:6);
SCTP=quantile(WCovX(:,2),0.82); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxBAHR == 1); idxLR=find(idxBALR == 1); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2)); YH=fH; XH=xH; 
XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']);
disp(['(1.3*) Controversy of BACH2''s role - GSE32918']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1); 
nowDthRt=0; bom=[]; Achk=[];
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
HRcolor=HRcolor2; HRcolorLine=HRcolorLine2;  
LRcolor=LRcolor2; LRcolorLine=LRcolorLine2; 
fsTt=10; LW = 2; PN1=plot(-0.5,-0.5,'color','w'); hold on;
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG1=plot(spHRt,spHRY,':','LineWidth',LW,'Color',HRcolorLine);
LRG1=plot(spLRt,spLRY,':','LineWidth',LW,'Color',LRcolorLine);


%%%whole GSE10846; for BACH2
%%%
Tobs=sDLBCL4(:,1); Delta=sDLBCL4(:,2); Icensor=sDLBCL4(:,3); WCovX=sDLBCL4(:,4:6);
SCTP=quantile(WCovX(:,2),0.85); idxBALR=[ (WCovX(:,2) > SCTP ) ]*1; idxBAHR=[ (WCovX(:,2) <= SCTP ) ]*1;
idxHR=find(idxBAHR == 1); idxLR=find(idxBALR == 1); HRd=Delta(idxHR ); HRt=Tobs(idxHR ); LRd=Delta(idxLR ); LRt=Tobs(idxLR ); 
AHH = [HRt.*((HRt <= 60)+0)+60*((HRt > 60)+0) 1-HRd.*((HRt <= 60)+0)]; BHH = sortrows(AHH,1);
ALL = [LRt.*((LRt <= 60)+0)+60*((LRt > 60)+0) 1-LRd.*((LRt <= 60)+0)]; BLL = sortrows(ALL,1);
[fH,xH] = ecdf(BHH(:,1),'censoring',BHH(:,2)); [fL,xL] = ecdf(BLL(:,1),'censoring',BLL(:,2)); YH=fH; XH=xH; 
XH=[xH; 60]; NN=size(fH); YH=[fH; fH(NN(1))]; XL=xL; YL=fL; XL=[xL; 60]; NN=size(fL); YL=[fL; fL(NN(1))]; disp(['   ']);
disp(['(1.4*) Controversy of BACH2''s role - GSE10846']); pvLR=logrankLF(BHH', BLL'); disp(['   ']); disp(['   ']);

%figures
csHRt=sort(HRt( find(HRd == 0) )); csHRY=0*csHRt+1; nnn=sum(XH > -1);
nowDthRt=0; bom=[]; Achk=[];
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
HRcolor=HRcolor1; HRcolorLine=HRcolorLine1;  
LRcolor=LRcolor1; LRcolorLine=LRcolorLine1; fsTt=10; LW = 2; 
PN2=plot(-0.6,-0.6, '.', 'color', 'w'); PNc=plot(-0.55,-0.55, '.', 'color', 'w');
plot(csHRt,csHRY,'+','Color',HRcolor); plot(csLRt,csLRY,'+','Color',LRcolor);
HRG2=plot(spHRt,spHRY,'LineWidth',LW,'Color',HRcolorLine); LRG2=plot(spLRt,spLRY,'LineWidth',LW,'Color',LRcolorLine); hold off;
LV=legend([HRG2 LRG2 PNc PN2 HRG1 LRG1], {'BACH2^{-} (n = 350)', 'BACH2^{+} (n = 62)',  '     ', '     ',  'BACH2^{-} (n = 137)', 'BACH2^{+} (n = 30)'}, 'Location', 'southwest'); LV=legend('boxoff'); set(LV,'FontSize',11);
text(2.5,0.35, 'GSE10846 (p = 0.0253)','FontWeight','bold','FontSize',11,'HorizontalAlignment','left'); 
text(2.5,0.16, 'GSE32918 (p = 0.303)','FontWeight','bold','FontSize',11,'HorizontalAlignment','left'); 
set(gcf,'Color',[1,1,1]); axis tight; ylim([0 1.05]); xlim([0 60]);
ylabel(['Survival  Rate']); xlabel(['Months']);
title(['GSE32918 - GSE10846'],'FontSize',fsTt, 'FontWeight', 'bold');  



