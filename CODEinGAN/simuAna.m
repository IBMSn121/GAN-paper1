function simuAna(c2ds,chks)
tc2d=c2ds; check=chks; InitialTime=datevec(now); CTP0=datevec(now); Cnum=1; Ctimes=10*[1:100]; JA=[400 800 1600 3200]; 
ngg=size(JA); ngg=ngg(2); tJ=max(JA); p=10; q=5; p0=0; Vrrr=[1.5]; Vc2d=[tc2d]; ddd=2; emso=0.2; pi1=0.4; pi2=0.5; pi=[pi1 pi2]; Rpt1=[]; Rpt2=[]; Rpt3=[]; disp(['Starting Simulations at   ' num2str(InitialTime(4)) ' : ' num2str(InitialTime(5))] ); for gggg = 1 : 1; rrr=Vrrr(gggg);
for g = 1 : 1; c2d=Vc2d(g); A1=[]; A2=[]; A3=[]; A4=[]; Apd1=[]; Apd2=[]; Apd3=[]; Apd4=[]; Afdr1=[]; Afdr2=[]; Afdr3=[]; Afdr4=[];
for ggg = 1 : check
BBn=( binornd(1,pi1,p-1,q).*(2*binornd(1,pi2,p-1,q)-1) ).*unifrnd(0.5,1,p-1,q);
BB=[BBn; 0.9999*ones(1,q)];
TWp1=chi2rnd(ddd,tJ,p-1); TW1=chi2rnd(ddd,tJ,1); TWp1=TWp1+rrr*TW1*ones(1,p-1); TW=[TWp1 TW1];
bme=chi2rnd(c2d,tJ,p); res=normrnd(0,sqrt(2*c2d),tJ,q); TT=chi2rnd(c2d,tJ,q); tT=chi2rnd(c2d,tJ,p);
X=TW+bme; TtV=(TT).^2 + 2*sin(TT); TY0=TtV+(X+cos(tT))*BB+res; for gg = 1 : ngg
J=JA(gg); W=TW(1:J,:); Y0=TY0(1:J,:); BpLSE=inv(W'*W)*W'*Y0; 
BpLaso=BpLSE; BpRR=BpLSE; BAWTE=BpLSE; for tt=1:q;
BpLaso(:,tt)=regLasso1(W, Y0(:,tt)); BpRR(:,tt)=rrHKB(W, Y0(:,tt));
[Beta,Pv] = regAWTE(W, Y0(:,tt)); BAWTE(:,tt)=Beta; end; AWTE=BAWTE; 
TE=(abs(BB)>0); veee=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45]; aaww=ones(1,9); for jj=1:9; emso=veee(jj); 
IEAWTE=(abs(AWTE)>emso); IEtAWTE=(IEAWTE>0&TE>0); NEtAWTE=(IEAWTE==0&TE==0);
aaww(jj)=(sum( sum(IEtAWTE) )/sum( sum(TE) )) + (sum( sum(NEtAWTE) )/sum( sum(1-TE) )); end; emso=veee(find(aaww==max(aaww),1)); 
IEAWTE=(abs(AWTE)>emso); IEtAWTE=(IEAWTE>0&TE>0);
%%%
aaww=ones(1,9); for jj=1:9; emso=veee(jj); IEBpLSE=(abs(BpLSE)>emso); IEtBpLSE=(IEBpLSE>0&TE>0); NEtBpLSE=(IEBpLSE==0&TE==0);
aaww(jj)=(sum( sum(IEtBpLSE) )/sum( sum(TE) )) + (sum( sum(NEtBpLSE) )/sum( sum(1-TE) )); end; emso=veee(find(aaww==max(aaww),1));
IEBpLSE=(abs(BpLSE)>emso); IEtBpLSE=(IEBpLSE>0&TE>0); 
%%%
aaww=ones(1,9); for jj=1:9; emso=veee(jj); IEBpRR=(abs(BpRR)>emso); IEtBpRR=(IEBpRR>0&TE>0); NEtBpRR=(IEBpRR==0&TE==0);
aaww(jj)=(sum( sum(IEtBpRR) )/sum( sum(TE) )) + (sum( sum(NEtBpRR) )/sum( sum(1-TE) )); end; emso=veee(find(aaww==max(aaww),1));
IEBpRR=(abs(BpRR)>emso); IEtBpRR=(IEBpRR>0&TE>0);
%%%
aaww=ones(1,9); for jj=1:9; emso=veee(jj); IEBpLaso=(abs(BpLaso)>emso); IEtBpLaso=(IEBpLaso>0&TE>0); NEtBpLaso=(IEBpLaso==0&TE==0);
aaww(jj)=(sum( sum(IEtBpLaso) )/sum( sum(TE) )) + (sum( sum(NEtBpLaso) )/sum( sum(1-TE) )); end; emso=veee(find(aaww==max(aaww),1));
IEBpLaso=(abs(BpLaso)>emso); IEtBpLaso=(IEBpLaso>0&TE>0);
%%%
aaa=[sum( sum(abs(BpLSE-BB)) ) sum( sum(abs(BpRR-BB)) ) sum( sum(abs(BB-AWTE)) ) sum( sum(abs(BB-BpLaso)) )];
%PD
aPD=[sum( sum(IEtBpLSE) )/sum( sum(TE) ) sum( sum(IEtBpRR) )/sum( sum(TE) ) sum( sum(IEtAWTE) )/sum( sum(TE) ) sum( sum(IEtBpLaso) )/sum( sum(TE) )];
%FDR
aFDR=1-[sum( sum(IEtBpLSE) )/sum( sum(IEBpLSE) ) sum( sum(IEtBpRR) )/sum( sum(IEBpRR) ) sum( sum(IEtAWTE) )/sum( sum(IEAWTE) ) sum( sum(IEtBpLaso) )/sum( sum(IEBpLaso) )];
if gg == 1
A1=[A1; aaa ];
Apd1=[Apd1; aPD];
Afdr1=[Afdr1; aFDR];
end
if gg == 2
A2=[A2; aaa];
Apd2=[Apd2; aPD];
Afdr2=[Afdr2; aFDR];
end
if gg == 3
A3=[A3; aaa];
Apd3=[Apd3; aPD];
Afdr3=[Afdr3; aFDR];
end  
if gg == 4
A4=[A4; aaa];
Apd4=[Apd4; aPD];
Afdr4=[Afdr4; aFDR];
end 
end %gg

if ggg > Ctimes(Cnum)
Cnum=Cnum+1; CTP1=datevec(now); CTP=(CTP1-CTP0)*[0 0 86400 3600 60 1]'/60; UTS='mins';
disp(['     Using' ': ' num2str(CTP) ' ' UTS ' ; current run = ' num2str(ggg) ' /' num2str(check) ])
end
end %ggg

if gggg == 1
Rpt1=[Rpt1; mean(A1); mean(A2); mean(A3); mean(A4)];
Rpt2=[Rpt2; mean(Apd1); mean(Apd2); mean(Apd3); mean(Apd4)]; Rpt3=[Rpt3; mean(Afdr1); mean(Afdr2); mean(Afdr3); mean(Afdr4)];
end
end %g
CTP1=datevec(now); CTP=(CTP1-CTP0)*[0 0 86400 3600 60 1]'/60; UTS='mins';
disp(['   Totally Using' ': ' num2str(CTP) ' ' UTS ' with Sigma2 = ' num2str(tc2d) ' and total runs = ' num2str(check)])
end %gggg


%%%PLOT
AAA=([A1 A2 A3 A4]); 
Name_S=char('\sigma^2 = ', '      r = ', '     \pi = '); scenario_num=1; Vvvv=Vc2d; Ftsz=10;
%Histogram for INER
ck=48; cc=[0:1/ck:1]'; cc=[cc cc cc]; cclr2=cc(41,:); cclr1=cclr2;
cc=hsv(ck); cclr0=cc(44,:);  cclr1(2)=0.65; cclr3=cclr2; cclr3(1)=0.75; ttt=0.2;

ccd=4; bw=2.5; uboX=100;
for tt=1:1
for gg = 1 : 4
subplot(ccd,4,gg+4*tt-4)
c = colormap(lines(4));
ttt=gg; h = histogram(AAA(:,3+ccd*(gg-1))); td=0.8+ttt; tc=0.94; h.FaceColor = c(3,:); h.EdgeColor = [tc tc tc]; h.BinWidth = bw;
xlim([0 uboX]); hold on; set(gcf,'Color',[1,1,1]);

x = AAA(:,4+ccd*(gg-1)); %LASSO
h1 = histogram(x(find(x<uboX))); td=0.6+ttt; tc=0.94;
h1.FaceColor = c(4,:);
h1.EdgeColor = [tc tc tc];
h1.BinWidth = bw;
xlim([0 uboX]);

x = AAA(:,2+ccd*(gg-1)); %RRE
h1 = histogram(x(find(x<uboX))); td=0.4+ttt; tc=0.94;
h1.FaceColor = c(2,:);
h1.EdgeColor = [tc tc tc];
h1.BinWidth = bw;
xlim([0 uboX]);

x = AAA(:,1+ccd*(gg-1)); %LSE
h1 = histogram(x(find(x<uboX))); td=0.2+ttt; tc=0.94;
h1.FaceColor = c(1,:);
h1.EdgeColor = [tc tc tc];
h1.BinWidth = bw;
xlim([0 uboX]);

if gg == 1; ylabel(['Histogram of INER'], 'FontWeight','bold','FontName','Time New Roman','FontSize',8); title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 2; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 3; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 4; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;

end
end
legend({'AWTE', 'LASSO', 'RRE', 'LSE' });
hold off;


%Scatter plot for FDR = X, PD = Y;
APD=([Apd1 Apd2 Apd3 Apd4]); 
AFDR=([Afdr1 Afdr2 Afdr3 Afdr4]);
ccd=4; c = colormap(lines(4)); lboY=0.5; uboX=0.9;
for tt=2:2
for gg = 1 : 4
subplot(ccd,4,gg+4*tt-4)
scatter( AFDR(:,3+ccd*(gg-1)), APD(:,3+ccd*(gg-1)), 'Marker', 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', c(3,:), 'SizeData', 5  );
%AWTE
hold on
set(gcf,'Color',[1,1,1]);
xlim([0 uboX])
ylim([lboY 1])
scatter( AFDR(:,4+ccd*(gg-1)), APD(:,4+ccd*(gg-1)), 'Marker', '+', 'MarkerFaceColor', c(4,:), 'MarkerEdgeColor', c(4,:), 'SizeData', 5 );
xlim([0 uboX])
ylim([lboY 1])
%LASSO
scatter( AFDR(:,2+ccd*(gg-1)), APD(:,2+ccd*(gg-1)), 'Marker', '^', 'MarkerFaceColor', c(2,:), 'MarkerEdgeColor', c(2,:), 'SizeData', 2 );
xlim([0 uboX])
ylim([lboY 1])
%RRE
scatter( AFDR(:,1+ccd*(gg-1)), APD(:,1+ccd*(gg-1)), 'Marker', 's', 'MarkerFaceColor', c(1,:), 'MarkerEdgeColor', c(1,:), 'SizeData', 2 );
xlim([0 uboX])
ylim([lboY 1])
%LSE
if gg == 1; ylabel(['FDR in X, PD in Y'], 'FontWeight','bold','FontName','Time New Roman','FontSize',8); title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 2; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 3; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 4; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
end
end
legend({'AWTE', 'LASSO', 'RRE', 'LSE' });
hold off;


%Square root
yrd=1; lboY=-0.1; scenario_num=1; dimplt=5; ck=48; cc=[0:1/ck:1]'; cc=[cc cc cc]; cclr2=cc(41,:); cclr1=cclr2;
cc=hsv(ck); cclr0=cc(44,:);  cclr1(2)=0.65; cclr3=cclr2; cclr3(1)=0.75;
c = colormap(lines(4)); C=[ones(4,1)*c(1,:); ones(4,1)*c(2,:); ones(4,1)*c(3,:); ones(4,1)*c(4,:)];
ccd=4; bw=2.5; tt=3; Vvvv=Vc2d;
for gg = 1 : 4
subplot(ccd,4,gg+4*tt-4)
aaa=AFDR(:,3+ccd*(gg-1)); bbb=APD(:,3+ccd*(gg-1));
AAWTE=sqrt(  ( aaa ).^2 + ( 1-bbb ).^2  );

aaa=AFDR(:,4+ccd*(gg-1)); bbb=APD(:,4+ccd*(gg-1));
AALaso=sqrt(  ( aaa ).^2 + ( 1-bbb ).^2  );

aaa=AFDR(:,2+ccd*(gg-1)); bbb=APD(:,2+ccd*(gg-1));
AARR=sqrt(  ( aaa ).^2 + ( 1-bbb ).^2  );

aaa=AFDR(:,1+ccd*(gg-1)); bbb=APD(:,1+ccd*(gg-1));
AALSE=sqrt(  ( aaa ).^2 + ( 1-bbb ).^2  );

c0=zeros(4,3);
if tt == 3; boxplot([AAWTE AALaso AARR AALSE] ,'colors', c0, 'labels',{' AW',' LA',' RR',' LS'}, 'OutlierSize',6, 'symbol', '.'); end;
hold on
set(gcf,'Color',[1,1,1]);
ylim([lboY yrd])
if gg == 1; ylabel(['Square root of cFP'], 'FontWeight','bold','FontName','Time New Roman','FontSize',8,'HorizontalAlignment','center'); title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 2; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 3; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
if gg == 4; title(['sample size = ' num2str(JA(gg))], 'FontWeight','bold','FontName','Time New Roman','FontSize',Ftsz); end;
end

tt=4;
subplot(ccd,4,3+4*tt-4)
title(['*In the case of ' Name_S(scenario_num,:) num2str(Vvvv(1)) ', where cFP=FDR^2+(1-PD)^2'], 'FontWeight','bold','FontName','Time New Roman','FontSize',11,'Color',[0.749019622802734 0 0.749019622802734],'HorizontalAlignment','right'); 
set(gca,'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
hold off;





