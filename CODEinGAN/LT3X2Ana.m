function LT3X2Ana(dX1,dY1,dX2,dY2,dX3,dY3,dX4,dY4,dX60,dY60)

%%%DLBCL GSE4475
%BCL6	SPIB	BACH2	IRF8	OCT2
DTrio=dX1;
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
DTO=dY1;
ganX=5; ganY=10;
%%%%%%%%%%%%
%Whole group mimus;
Trendlse=zeros(ganX,ganY); TrendlseB=zeros(ganX,ganY); TrendlseU=zeros(ganX,ganY); TrendlseL=zeros(ganX,ganY); TrendP=zeros(ganX,ganY);
for ggg = 1 : ganX
for gg = 1 : ganY

Xg=DTrio(:,ggg); Yg=DTO(:,gg); Tobs=abs(Yg)+1; 
Wr=Xg; Y0=Yg; J=sum((Tobs > 0));
Is=eye(J)-ones(J)/J; bpLSEfct=inv(Wr'*Is*Wr)*Wr'*Is*Y0;
Trendlse(ggg,gg)=bpLSEfct; 
[bpLSE1,bint,r,rint,stats] = regress(Y0,[Wr*0+1 Wr]); TrendlseB(ggg,gg)=bpLSE1(2); TrendlseU(ggg,gg)=bint(2,2); TrendlseL(ggg,gg)=bint(2,1);
TrendP(ggg,gg)=stats(3);
end
end
AAA=[Trendlse(1,[4 6 9])'; Trendlse(2,[3 6 8])'; Trendlse(3,[5:10])'; Trendlse(4,[5:7])'; Trendlse(5,[4 6 8])'];
Trendlsea=TrendlseB;
AAAB=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])']; 

Trendlsea=TrendlseU;
AAAU=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendlseL;
AAAL=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendP;
AAAP=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

AAA1P=AAAP;
AAA1B=AAAB; AAA1U=AAAU; AAA1L=AAAL;
%


%%%DLBCL GSE32918
%BCL6	SPIB	BACH2	IRF8	OCT2
DTrio=dX2;
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
DTO=dY2;
ganX=5; ganY=10;
%%%%%%%%%%%%
%Whole group mimus;
Trendlse=zeros(ganX,ganY); TrendlseB=zeros(ganX,ganY); TrendlseU=zeros(ganX,ganY); TrendlseL=zeros(ganX,ganY); TrendP=zeros(ganX,ganY);
for ggg = 1 : ganX
for gg = 1 : ganY

Xg=DTrio(:,ggg); Yg=DTO(:,gg); Tobs=abs(Yg)+1;   
Wr=Xg; Y0=Yg; J=sum((Tobs > 0));
Is=eye(J)-ones(J)/J; bpLSEfct=inv(Wr'*Is*Wr)*Wr'*Is*Y0;
Trendlse(ggg,gg)=bpLSEfct; 
[bpLSE1,bint,r,rint,stats] = regress(Y0,[Wr*0+1 Wr]); TrendlseB(ggg,gg)=bpLSE1(2); TrendlseU(ggg,gg)=bint(2,2); TrendlseL(ggg,gg)=bint(2,1);
TrendP(ggg,gg)=stats(3);
end
end
AAA=[Trendlse(1,[4 6 9])'; Trendlse(2,[3 6 8])'; Trendlse(3,[5:10])'; Trendlse(4,[5:7])'; Trendlse(5,[4 6 8])'];
Trendlsea=TrendlseB;
AAAB=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])']; 

Trendlsea=TrendlseU;
AAAU=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendlseL;
AAAL=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendP;
AAAP=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

AAA2P=AAAP;
AAA2B=AAAB; AAA2U=AAAU; AAA2L=AAAL; 


%%%DLBCL GSE31312
%BCL6	SPIB	BACH2	IRF8	OCT2
DTrio=dX3;
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
DTO=dY3;
ganX=5; ganY=10;
%%%%%%%%%%%%
%Whole group mimus;
Trendlse=zeros(ganX,ganY); TrendlseB=zeros(ganX,ganY); TrendlseU=zeros(ganX,ganY); TrendlseL=zeros(ganX,ganY); TrendP=zeros(ganX,ganY);
for ggg = 1 : ganX
for gg = 1 : ganY

Xg=DTrio(:,ggg); Yg=DTO(:,gg); Tobs=abs(Yg)+1;   
Wr=Xg; Y0=Yg; J=sum((Tobs > 0));
Is=eye(J)-ones(J)/J; bpLSEfct=inv(Wr'*Is*Wr)*Wr'*Is*Y0;
Trendlse(ggg,gg)=bpLSEfct; 
[bpLSE1,bint,r,rint,stats] = regress(Y0,[Wr*0+1 Wr]); TrendlseB(ggg,gg)=bpLSE1(2); TrendlseU(ggg,gg)=bint(2,2); TrendlseL(ggg,gg)=bint(2,1);
TrendP(ggg,gg)=stats(3);
end
end
AAA=[Trendlse(1,[4 6 9])'; Trendlse(2,[3 6 8])'; Trendlse(3,[5:10])'; Trendlse(4,[5:7])'; Trendlse(5,[4 6 8])'];
Trendlsea=TrendlseB;
AAAB=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])']; 

Trendlsea=TrendlseU;
AAAU=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendlseL;
AAAL=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendP;
AAAP=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

AAA3P=AAAP;
AAA3B=AAAB; AAA3U=AAAU; AAA3L=AAAL; 
%


%%%DLBCL GSE10846
%BCL6	SPIB	BACH2	IRF8	OCT2
DTrio=dX4;
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
DTO=dY4;
ganX=5; ganY=10;
%%%%%%%%%%%%
%Whole group mimus;
Trendlse=zeros(ganX,ganY); TrendlseB=zeros(ganX,ganY); TrendlseU=zeros(ganX,ganY); TrendlseL=zeros(ganX,ganY); TrendP=zeros(ganX,ganY);
for ggg = 1 : ganX
for gg = 1 : ganY

Xg=DTrio(:,ggg); Yg=DTO(:,gg); Tobs=abs(Yg)+1;   
Wr=Xg; Y0=Yg; J=sum((Tobs > 0));
Is=eye(J)-ones(J)/J; bpLSEfct=inv(Wr'*Is*Wr)*Wr'*Is*Y0;
Trendlse(ggg,gg)=bpLSEfct; 
[bpLSE1,bint,r,rint,stats] = regress(Y0,[Wr*0+1 Wr]); TrendlseB(ggg,gg)=bpLSE1(2); TrendlseU(ggg,gg)=bint(2,2); TrendlseL(ggg,gg)=bint(2,1);
TrendP(ggg,gg)=stats(3);
end
end
AAA=[Trendlse(1,[4 6 9])'; Trendlse(2,[3 6 8])'; Trendlse(3,[5:10])'; Trendlse(4,[5:7])'; Trendlse(5,[4 6 8])'];
Trendlsea=TrendlseB;
AAAB=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])']; 

Trendlsea=TrendlseU;
AAAU=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendlseL;
AAAL=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendP;
AAAP=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

AAA4P=AAAP;
AAA4B=AAAB; AAA4U=AAAU; AAA4L=AAAL;


%GSE60_tumor only
%BCL6	SPIB	BACH2	IRF8	OCT2
DTrio=dX60;
%p21	MYC	p53	BCL2	NFKB1	IRF4	Blimp1	AID	p27	ATR
DTO=dY60;
ganX=5; ganY=10;
%%%%%%%%%%%%
%Whole group mimus;
Trendlse=zeros(ganX,ganY); TrendlseB=zeros(ganX,ganY); TrendlseU=zeros(ganX,ganY); TrendlseL=zeros(ganX,ganY); TrendP=zeros(ganX,ganY);
for ggg = 1 : ganX
for gg = 1 : ganY

Xg=DTrio(:,ggg); Yg=DTO(:,gg); Tobs=abs(Yg)+1;   
Wr=Xg; Y0=Yg; J=sum((Tobs > 0));
Is=eye(J)-ones(J)/J; bpLSEfct=inv(Wr'*Is*Wr)*Wr'*Is*Y0;
Trendlse(ggg,gg)=bpLSEfct; 
[bpLSE1,bint,r,rint,stats] = regress(Y0,[Wr*0+1 Wr]); TrendlseB(ggg,gg)=bpLSE1(2); TrendlseU(ggg,gg)=bint(2,2); TrendlseL(ggg,gg)=bint(2,1);
TrendP(ggg,gg)=stats(3);
end
end
AAA=[Trendlse(1,[4 6 9])'; Trendlse(2,[3 6 8])'; Trendlse(3,[5:10])'; Trendlse(4,[5:7])'; Trendlse(5,[4 6 8])'];
Trendlsea=TrendlseB;
AAAB=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])']; 

Trendlsea=TrendlseU;
AAAU=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendlseL;
AAAL=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

Trendlsea=TrendP;
AAAP=[Trendlsea(1,[4 6 9])'; Trendlsea(2,[3 6 8])'; Trendlsea(3,[5:10])'; Trendlsea(4,[5:7])'; Trendlsea(5,[4 6 8])'];

AAA60P=AAAP; AAA60B=AAAB; AAA60U=AAAU; AAA60L=AAAL;

AAAP=[AAA60P AAA1P AAA2P AAA3P AAA4P];
AAAB=[AAA60B AAA1B AAA2B AAA3B AAA4B]; AAAU=[AAA60U AAA1U AAA2U AAA3U AAA4U]; AAAL=[AAA60L AAA1L AAA2L AAA3L AAA4L];


%plot draw here; only main6
UpRC0lor=[0.95686274766922 0.282352954149246 0.61176472902298]; DnRC0lor=[0.0352941192686558 0.878431379795074 0.670588254928589];
GC0lor=[0.494117647409439 0.494117647409439 0.494117647409439];
GC2lor=[0 0.627451002597809 0.988235294818878]; GC2lora=[0.831372559070587 0.901960790157318 0.945098042488098];
tGC2lor=[0.972549021244049 0.592156887054443 0.10196078568697];
tGC2lorC=[0.972549021244049 0.690196096897125 0.200000002980232];
tGC2lora=[0.87058824300766 0.760784327983856 0.615686297416687];
ttGC2lor=[1 0.39215686917305 0.156862750649452]; 
ttGC2lora=[0.972549021244049 0.796078443527222 0.725490212440491]; 
GC3lor=[0.749019622802734 0 0.513725519180298]; GC3lora=[0.607843160629272 0.447058826684952 0.556862771511078];
GC3lorb=[0.658823549747467 0.5686274766922 0.627451002597809];
GC3lorC=[0.749019622802734 0.321568638086319 0.650980412960052];
GGreen1=[0.898039221763611 1 0.933333337306976]; 
GGreen4=[0 0.800000011920929 0.266666680574417];
GPurple1=[0.988235294818878 0.921568632125854 0.933333337306976];
GPurple4=[0.874509811401367 0.121568627655506 0.243137255311012];
%
txmain6={'SPIB-IRF4', 'SPIB-AID','BACH2-IRF4','BACH2-AID','OCT2-IRF4', 'OCT2-AID'}; idmain6=[5 6 8 10 17 18]; 
idxVAR=idmain6;

Xempty=3; Yempty=1; xLBnd=0; xRBnd=22; fntsz=9; yRBnd1=-5; yRBnd2=15;   slightA=0.1; xmksz=10; oclr=0.4;
plot(xRBnd,6+Yempty,'.','Color',[1 1 1],'MarkerSize',80); hold on;
EptFix=[0 0 -0.5 -0.5 0 0.5]; EptRng=[-0.5 1.5];

USgeneS={'SPIB','BACH2','OCT2'}; DSgeneS={'IRF4','AID'};
Xmain=1; Ymain=1; RRmk=[1 0.8 0.6 0.4 0.2]; fntsz=8;  
CenterX=0*[1:6];
ccum=0; JmpInt=2.5; textCor=[1 1 1]*0.933333337306976; StatX=1;
for ttt=1:3;
for ggg = 1 :2;
ypstn=10; rngX=[StatX 1.25*StatX+JmpInt]+ccum*JmpInt;
xmksz=7+Xmain*15; xlnsz=2.5+Xmain*0.5; mecR=0.5-Xmain*0.2; Afntsz=fntsz+Xmain*2; 
adjY=1.6; adjX=-0.38; if ttt == 2; adjX=-0.39; end;

for tts = 1 : 4
tt=5-tts;
plot(mean(rngX)+adjX+tt*slightA/4,ypstn+adjY-tt*slightA/4,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]*(0.5+0.1*tt),'MarkerSize',xmksz,'Marker','diamond','LineWidth',xlnsz);
end
plot(mean(rngX)+adjX,ypstn+adjY,'MarkerFaceColor',tGC2lorC,'MarkerEdgeColor',tGC2lorC,'MarkerSize',xmksz,'Marker','diamond','LineWidth',xlnsz);
text(mean(rngX)+adjX-0.04,ypstn+adjY, USgeneS(ttt),'FontWeight','bold','FontSize',Afntsz,'HorizontalAlignment','center','Color',[1 1 1]*0.25); 
xmksz=10+Ymain*15; xlnsz=0.5+Ymain*0.5; mecR=0.5-Ymain*0.2; Afntsz=fntsz+Ymain*2;
adjY=1.6; adjX=0.38; if ttt == 2; adjX=0.39; end;
for tts = 1 : 4
tt=5-tts;
plot(mean(rngX)+adjX+tt*slightA/4,ypstn+adjY-tt*slightA/4,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1]*(0.5+0.1*tt),'MarkerSize',xmksz,'Marker','o','LineWidth',xlnsz);
end
plot(mean(rngX)+adjX,ypstn+adjY, 'MarkerFaceColor',GC3lorC,'MarkerEdgeColor',GC3lorC,'MarkerSize',xmksz+2.5,'Marker','o','LineWidth',xlnsz);
text(mean(rngX)+adjX+0.02,ypstn+adjY, DSgeneS(ggg),'FontWeight','bold','FontSize',Afntsz,'HorizontalAlignment','center','Color',textCor); 
ccum=ccum+1; CenterX(ccum)=mean(rngX);
xxL=CenterX(ccum)-1.1:0.2:CenterX(ccum)+1.1; yyL=0*xxL+4; plot(xxL,yyL,'k','LineWidth',2); eRng=EptRng+EptFix(ccum);
yyP=4:0.1:4.4; xxP=xxL(1)+0*yyP+0.1; plot(xxP,yyP,'Color',[1 1 1]*0.25,'LineWidth',1.5); 
text(xxP(1),4-0.4, num2str(eRng(1)),'FontWeight','bold','FontSize',Afntsz,'HorizontalAlignment','center','Color',[1 1 1]*0.25);
ptzero=xxP(1)-eRng(1);
yyP=4:0.1:4.4; xxP=xxL(sum(xxL>0))+0*yyP-0.1; plot(xxP,yyP,'Color',[1 1 1]*0.25,'LineWidth',1.5); 
text(xxP(1),4-0.4, num2str(eRng(2)),'FontWeight','bold','FontSize',Afntsz,'HorizontalAlignment','center','Color',[1 1 1]*0.25);
Yy=4:0.5:9.5; Xx=Yy*0+ptzero; plot(Xx,Yy,':','Color',[0.662745118141174 0.184313729405403 0.556862771511078],'LineWidth',3);


Tchk=5;
cchk=sum(AAAB(idxVAR(ccum),:) > 0); Icchk=(cchk==Tchk)+(cchk==0); GCclr=GC0lor;
if cchk==Tchk; GCclr=UpRC0lor; end; if cchk==0; GCclr=DnRC0lor; end;

Xx=AAAB(idxVAR(ccum),:)+ptzero; Yy=[5:9]; 
if Icchk == 1; plot(Xx,Yy,'-.','Color',GCclr,'LineWidth',3); end;
tt=1; ciLth=1.25;
YyC=Yy(tt)+[-2:2]/10; 
XxU=YyC*0+AAAU(idxVAR(ccum),tt)+ptzero; XxL=YyC*0+AAAL(idxVAR(ccum),tt)+ptzero;
plot(XxU,YyC,'Color',0.31*[1 1 1],'LineWidth',ciLth); plot(XxL,YyC,'Color',0.31*[1 1 1],'LineWidth',ciLth);

XxC=[AAAL(idxVAR(ccum),tt):(AAAU(idxVAR(ccum),tt)-AAAL(idxVAR(ccum),tt))/5:AAAU(idxVAR(ccum),tt)]+ptzero; YyC=XxC*0+Yy(tt); 
plot(XxC,YyC,'Color',0.0*[1 1 1],'LineWidth',ciLth);

if AAAB(idxVAR(ccum),tt)>0; GCclr=UpRC0lor; end; if AAAB(idxVAR(ccum),tt)<0; GCclr=DnRC0lor; end;
plot(AAAB(idxVAR(ccum),tt)+ptzero,Yy(tt),'Color',GCclr,'Marker','*','MarkerSize',7,'LineWidth',1);

%4 CI
ciLth=1.25; 
for tt = 2 : 5;
YyC=Yy(tt)+[-2:2]/10; 
XxU=YyC*0+AAAU(idxVAR(ccum),tt)+ptzero; XxL=YyC*0+AAAL(idxVAR(ccum),tt)+ptzero;
plot(XxU,YyC,'Color',0.31*[1 1 1],'LineWidth',ciLth); plot(XxL,YyC,'Color',0.31*[1 1 1],'LineWidth',ciLth);

XxC=[AAAL(idxVAR(ccum),tt):(AAAU(idxVAR(ccum),tt)-AAAL(idxVAR(ccum),tt))/5:AAAU(idxVAR(ccum),tt)]+ptzero; YyC=XxC*0+Yy(tt); 
plot(XxC,YyC,'Color',0.0*[1 1 1],'LineWidth',ciLth);

mkSet='square'; if tt == 3;  mkSet='diamond'; elseif tt == 4;  mkSet='^'; elseif tt == 5;  mkSet='o'; end;

if AAAB(idxVAR(ccum),tt)>0; GCclr=UpRC0lor; end; if AAAB(idxVAR(ccum),tt)<0; GCclr=DnRC0lor; end;
plot(AAAB(idxVAR(ccum),tt)+ptzero,Yy(tt),'Color',GCclr,'Marker',mkSet,'MarkerSize',7,'LineWidth',1);
plot(AAAB(idxVAR(ccum),tt)+ptzero,Yy(tt),'Marker',mkSet,'MarkerSize',5,'LineWidth',1,'MarkerFaceColor',GCclr*0.77,'MarkerEdgeColor',[1 1 1]);
end;


end;%ggg 
end;%ttt

Xswch=16.5; lsy=11; 
GseDB={'GSE4475', 'GSE32918', 'GSE31312', 'GSE10846'}; 
Afntsz=11;  textClr=[0 0.447058826684952 0.74117648601532]; 
GCclr=GC0lor; sdd=1; 
for tt = 1 : 4
mkSet='square'; if tt == 2;  mkSet='diamond'; elseif tt == 3;  mkSet='^'; elseif tt == 4;  mkSet='o'; end;
plot(Xswch,lsy-tt*sdd,'Color',GCclr,'Marker',mkSet,'MarkerSize',7,'LineWidth',1);
plot(Xswch,lsy-tt*sdd,'Marker',mkSet,'MarkerSize',5,'LineWidth',1,'MarkerFaceColor',GCclr*0.77,'MarkerEdgeColor',[1 1 1]);
text(Xswch+0.2,lsy-tt*sdd,GseDB(tt),'FontWeight','bold','FontSize',Afntsz, 'Color',textClr,'HorizontalAlignment','left','LineStyle','none');
end

tt=5; GCclr=UpRC0lor; Xx=Xswch+[-1 0 1]/10; Yy=[1 1 1]*(lsy-tt*sdd); 
plot(Xx,Yy,'-.','Color',GCclr,'LineWidth',2);  
text(Xx(2)+0.2,Yy(2),'Consistent up-regulation','FontWeight','bold','FontSize',Afntsz, 'Color',textClr,'HorizontalAlignment','left','LineStyle','none');

tt=6; GCclr=DnRC0lor; Xx=Xswch+[-1 0 1]/10; Yy=[1 1 1]*(lsy-tt*sdd); 
plot(Xx,Yy,'-.','Color',GCclr,'LineWidth',2);  
text(Xx(2)+0.2,Yy(2),'Consistent down-regulation','FontWeight','bold','FontSize',Afntsz, 'Color',textClr,'HorizontalAlignment','left','LineStyle','none');

plot(Xswch,lsy,'Color',GC0lor,'Marker','*','MarkerSize',7,'LineWidth',1);
text(Xswch+0.2,lsy,'GSE60','FontWeight','bold','FontSize',Afntsz, 'Color',textClr,'HorizontalAlignment','left','LineStyle','none');

hold off; set(gcf,'Color',[1,1,1]); axis tight; xlim([xLBnd xRBnd]); ylim([yRBnd1 yRBnd2]); 
text(mean(CenterX),2,'Slope  of  linear  trend','FontWeight','bold','FontSize',Afntsz+1,'HorizontalAlignment','center','LineStyle','none');



