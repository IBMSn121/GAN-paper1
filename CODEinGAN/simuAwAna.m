function simuAwAna(chks)
QQQ=[5 10 100 1000]; ngg=20; 
InitialTime=datevec(now); CTP0=datevec(now); Ctimes=20*[1:50]; JA=500*[1:ngg]; 
tJ=max(JA); p=10; p0=0; Vrrr=[1.5]; Vc2d=[1]; check=chks; ddd=2; emso=0.4; pi1=0.4; pi2=0.5; pi=[pi1 pi2]; 
Rpt1=[]; Rpt2=[]; Rpt3=[]; Rpt4=[]; disp(['Starting Simulations at   ' num2str(InitialTime(4)) ' : ' num2str(InitialTime(5))] ); A1=[]; A2=[]; A3=[]; A4=[];Apd1=[]; Apd2=[]; Apd3=[]; Apd4=[]; 
Afdr1=[]; Afdr2=[]; Afdr3=[]; Afdr4=[]; for gggg = 1 : 4; q=QQQ(gggg); Cnum=1; for g = 1 : 1;
c2d=Vc2d(g); rrr=Vrrr(g); for ggg = 1 : check; BBn=( binornd(1,pi1,p-1,q).*(2*binornd(1,pi2,p-1,q)-1) ).*unifrnd(0.5,1,p-1,q);
BB=[BBn; 0.9999*ones(1,q)]; TWp1=chi2rnd(ddd,tJ,p-1); TW1=chi2rnd(ddd,tJ,1); TWp1=TWp1+rrr*TW1*ones(1,p-1); TW=[TWp1 TW1]; bme=chi2rnd(c2d,tJ,p); res=normrnd(0,sqrt(2*c2d),tJ,q); TT=chi2rnd(c2d,tJ,q); tT=chi2rnd(c2d,tJ,p); X=TW+bme; TtV=(TT).^2 + 2*sin(TT); TY0=TtV+(X+cos(tT))*BB+res; aaa=[]; aPD=[]; aFDR=[]; extTest=randperm(JA(1)); extTest=extTest(1:JA(1)/10);
for gg = 1 : ngg; J=JA(gg); W=TW(1:J,:); Y0=TY0(1:J,:); k0=0.5; H0=W(:,p); g0=( H0>median(H0) ); cwp=((g0-k0*ones(J,1))'*H0); a0h=(g0-k0*ones(J,1))'*W/cwp; a0h(p)=0; [Beta,Pv] = regAWTE(W(extTest,:), Y0(extTest,1));
Wa=W-H0*ones(1,p)*diag(a0h); GWa=(Wa>ones(J,1)*median(Wa))-0.5*ones(J,p); AWTE=(GWa'*Y0)./(diag(GWa'*Wa)*ones(1,q));
Kstar=((g0-k0*ones(J,1))'*Y0)/cwp; AWTE(p,:)=Kstar-a0h*AWTE; TE=(abs(BB)>0); veee=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45]; KKTT=[]; aaww=ones(1,9); 
for jj=1:9; emso=veee(jj); IEAWTE=(abs(AWTE)>emso); IEtAWTE=(IEAWTE>0&TE>0); NEtAWTE=(IEAWTE==0&TE==0);
aaww(jj)=(sum( sum(IEtAWTE) )/sum( sum(TE) )) + (sum( sum(NEtAWTE) )/sum( sum(1-TE) )); end;
emso=veee(find(aaww==max(aaww),1)); aaa=[aaa sum( sum(abs(BB-AWTE)) )]; IEAWTE=(abs(AWTE)>emso); IEtAWTE=(IEAWTE>0&TE>0); aPD=[aPD sum( sum(IEtAWTE) )/sum( sum(TE) )]; aFDR=[aFDR (1-sum( sum(IEtAWTE) )/sum( sum(IEAWTE) ))]; end %gg
if gggg == 1; A1=[A1; aaa]; Apd1=[Apd1; aPD]; Afdr1=[Afdr1; aFDR]; end;
if gggg == 2; A2=[A2; aaa]; Apd2=[Apd2; aPD]; Afdr2=[Afdr2; aFDR]; end;
if gggg == 3; A3=[A3; aaa]; Apd3=[Apd3; aPD]; Afdr3=[Afdr3; aFDR]; end;  
if gggg == 4; A4=[A4; aaa]; Apd4=[Apd4; aPD]; Afdr4=[Afdr4; aFDR]; end; 
if ggg > Ctimes(Cnum)
Cnum=Cnum+1; CTP1=datevec(now); CTP=(CTP1-CTP0)*[0 0 86400 3600 60 1]'/60; UTS='mins';
disp(['     Using' ': ' num2str(CTP) ' ' UTS ' ; current run = ' num2str(ggg) ' /' num2str(check) ' ( number of down-stream genes q = ' num2str(q) ' )' ])
end
end %ggg
end %g
CTP1=datevec(now); CTP=(CTP1-CTP0)*[0 0 86400 3600 60 1]'/60; UTS='mins';
disp(['   Totally Using' ': ' num2str(CTP) ' ' UTS ' with total runs = ' num2str(check) ' ( number of down-stream genes q = ' num2str(q) ' )'])
end %gggg


ck=48; cc=[0:1/ck:1]'; cc=[cc cc cc]; cclr2=cc(41,:); cclr1=cclr2; 
cc=hsv(ck); cclr0=cc(44,:);  cclr1(2)=0.65; cclr3=cclr2; cclr3(1)=0.75; c = zeros(4,3); JA=0.5*[1:ngg];
subplot(2,6,1);
plot(JA, mean(log(A1)), 'color', c(1,:))
hold on;
plot(JA, mean(log(A2)), '--', 'color', c(2,:))
plot(JA, mean(log(A3)), ':', 'color',c(3,:))
plot(JA, mean(log(A4)), '*', 'color',c(4,:),'MarkerSize',3)
hold off;
set(gcf,'Color',[1,1,1]);
xlabel('Sample Size x 1000');
ylabel('ln INER');
axis tight;

subplot(2,6,2);
plot(JA, mean(Apd1), 'color', c(1,:))
hold on;
plot(JA, mean(Apd2), '--', 'color', c(2,:))
plot(JA, mean(Apd3), ':', 'color',c(3,:))
plot(JA, mean(Apd4), '*', 'color',c(4,:),'MarkerSize',3)
hold off;
set(gcf,'Color',[1,1,1]);
xlabel('Sample Size x 1000');
ylabel('PD');
axis tight;

subplot(2,6,3);
plot(JA, mean(Afdr1), 'color', c(1,:))
hold on;
plot(JA, mean(Afdr2), '--', 'color', c(2,:))
plot(JA, mean(Afdr3), ':', 'color',c(3,:))
plot(JA, mean(Afdr4), '*', 'color',c(4,:),'MarkerSize',3)
hold off;
set(gcf,'Color',[1,1,1]);
xlabel('Sample Size x 1000');
ylabel('FDR');
axis tight;
legend({'q =  5', 'q = 10', 'q = 100', 'q = 1000'});

subplot(2,6,[4 5]+3);
ck=48; cc=[0:1/ck:1]'; cc=[cc cc cc]; cclr2=cc(41,:); cclr1=cclr2;
cc=hsv(ck); cclr0=cc(44,:);  cclr1(2)=0.65; cclr3=cclr2; cclr3(1)=0.75;
JA=500*[1:ngg]; c = zeros(4,3); 
C=[ones(4,1)*c(1,:); ones(4,1)*c(2,:); ones(4,1)*c(3,:); ones(4,1)*c(4,:)];
plot(JA, mean(A1)/(10*5), 'color', c(1,:))
hold on;
plot(JA, mean(A2)/(10*10), '--', 'color', c(2,:))
plot(JA, mean(A3)/(10*100), ':', 'color',c(3,:))
plot(JA, mean(A4)/(10*1000), '*', 'color',c(4,:),'MarkerSize',3)
hold off;
set(gcf,'Color',[1,1,1]);
xlabel('Sample Size');
ylabel('INER∕pq');
axis tight;
legend({'q =  5', 'q = 10', 'q = 100', 'q = 1000' });

