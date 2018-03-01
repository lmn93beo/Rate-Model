clear all
rng(5)
T=500; dt=.5; T=T/dt;
n=2.2; %lin/nonlin. 1=lin, >1=n for nonlin
maxopt=0; %0 for nonlin
k=.01; xn=2;
Ne=400*xn; Npv=50*xn; Nvip=25*xn; Nsom=25*xn;
Ntot=Ne+Npv+Nvip+Nsom; % # of cell groups
taus=[20 10 10 10]; %time constants, e i


conprobs=[.15, .3, .3, .3,     1;
    .3, .3, 0.3, .3,     1;
    .3, .3, 0.01, .3,      0;
    .3, 0.01, .3, .01,      0];


W4=[0.0171   -0.9559   -0.0450   -0.5121;
    0.8535   -0.9900   -0.0900   -0.3073;
    1.9800   -0.1843         0   -0.7339;
    1.2848         0   -0.1399         0]; %2.2379 1.3276 (L=1.2076)


NCmat=zeros(size(W4));
NCmat(:,1)=Ne; NCmat(:,2)=Npv; NCmat(:,3)=Nvip; NCmat(:,4)=Nsom;
gsyns=W4./(NCmat.*conprobs(:,1:4));


NeI=[1:Ne]; NpvI=[Ne+1:Ne+Npv]; NvipI=[Ne+Npv+1:Ne+Npv+Nvip];  NsomI=[Ne+Nvip+Npv+1:Ne+Npv+Nvip+Nsom];
NIs{1}=NeI; NIs{2}=NpvI; NIs{3}=NvipI; NIs{4}=NsomI; NIs{5}=Ntot+1;



for ni=1 %over networks
   
    
    %%
    W=zeros(Ntot, Ntot);
    
    dvwn=9; %variance of weight equals 1/dvwn of mean
    
    for prei=1:4 %make weight matrix
        for posti=1:4
            if gsyns(posti,prei)~=0
                MU = log(abs(gsyns(posti,prei))^2 / sqrt((abs(gsyns(posti,prei))/dvwn)+abs(gsyns(posti,prei))^2));
                SIGMA = sqrt(log((abs(gsyns(posti,prei))/dvwn)/(abs(gsyns(posti,prei)))^2 + 1));
                
                rweights=lognrnd(MU, SIGMA, length(NIs{posti}), length(NIs{prei})).*sign(gsyns(posti,prei));
                if prei==1
                    rweights(rweights<0)=gsyns(posti,prei)/10; %clipping
                else
                    rweights(rweights>0)=gsyns(posti,prei)/10; %clipping
                end
                
                W(NIs{posti}, NIs{prei}) = (rand(length(NIs{posti}), length(NIs{prei})) < conprobs(posti, prei)).* rweights;
            end
            
        end
    end
    
    %%
    
    SpontRates=[12.8 24 8 8.9];
    
    InpRates=[93 74 0 0];
    
    BinpE=ones(1,T).*SpontRates(1); BinpPV=ones(1, T).*SpontRates(2); BinpVIP=ones(1, T).*SpontRates(3);  BinpSOM=ones(1, T).*SpontRates(4); %no reason for this to be over time....
    
    endT=700; stT=500;
    InpEk=(rand(Ne,1)>conprobs(1,end)); InpPVk=(rand(Npv,1)>conprobs(2,end)); InpVIPk=(rand(Nvip,1)>conprobs(3,end)); InpSOMk=(rand(Nsom,1)>conprobs(4,end)); %which cells dont get inp
    
    %run it
    DV=.5;
    Drates=[0 1 1 1]*5;
    vn=sqrt((Drates./DV)./.083);
    DepolsE=randn(Ne,1).*sqrt(Drates(1)/DV)+Drates(1); DepolsPV=randn(Npv,1).*sqrt(Drates(2)/DV)+Drates(2); DepolsVIP=randn(Nvip,1).*sqrt(Drates(3)/DV)+Drates(3); DepolsSOM=randn(Nsom,1).*sqrt(Drates(4)/DV)+Drates(4); DepolsE(DepolsE<0)=0; DepolsPV(DepolsPV<0)=0; DepolsVIP(DepolsVIP<0)=0; DepolsSOM(DepolsSOM<0)=0;
    
    
    for ti=1:10 %trials
        for aii=1:2 %ach or not
            InpE=zeros(Ne,T); InpPV=zeros(Npv,T); InpVIP=zeros(Nvip,T); InpSOM=zeros(Nsom,T);
            InpE(:,stT:endT)=InpE(:,stT:endT)+InpRates(1)+repmat(randn(Ne,1).*sqrt(2),1,endT-stT+1); %+repmat(randn(Ne,1).*InpRates(1)/40,1,endT-500+1);
            InpPV(:,stT:endT)=InpPV(:,stT:endT)+InpRates(2)+repmat(randn(Npv,1).*sqrt(2),1,endT-stT+1); %+repmat(randn(Npv,1).*InpRates(2)/40,1,endT-500+1);
            InpVIP(:,stT:endT)=InpVIP(:,stT:endT)+InpRates(3); %+repmat(randn(Nvip,1).*sqrt(InpRates(3)/40),1,endT-500+1);
            InpSOM(:,stT:endT)=InpSOM(:,stT:endT)+InpRates(4); %+repmat(randn(Nsom,1).*sqrt(InpRates(4)/40),1,endT-500+1);
            InpE(InpEk,:)=0; InpPV(InpPVk,:)=0; InpVIP(InpVIPk,:)=0; InpSOM(InpSOMk,:)=0;
            
            Re(:,1)=zeros(Ne,1); Rpv(:,1)=zeros(Npv,1); Rvip(:,1)=zeros(Nvip,1); Rsom(:,1)=zeros(Nsom,1);
            for t=2:T %over time
                dRe=(-Re(:,t-1)+k.*(max((InpE(:,t-1)+BinpE(1,t-1)+ DepolsE*(aii==1)+ W(NeI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(1);
                dRpv=(-Rpv(:,t-1)+k.*(max((InpPV(:,t-1)+BinpPV(1,t-1)+ DepolsPV*(aii==1)+W(NpvI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(2);
                dRvip=(-Rvip(:,t-1)+k.*(max((InpVIP(:,t-1)+BinpVIP(1,t-1)+ DepolsVIP*(aii==1)+ W(NvipI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(3);
                dRsom=(-Rsom(:,t-1)+k.*(max((InpSOM(:,t-1)+BinpSOM(1,t-1)+ DepolsSOM*(aii==1)+ W(NsomI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(4);
                
                
                Re(:,t)= Re(:,t-1)+dRe*dt;
                Rpv(:,t)= Rpv(:,t-1)+dRpv*dt;
                Rvip(:,t)= Rvip(:,t-1)+dRvip*dt;
                Rsom(:,t)= Rsom(:,t-1)+dRsom*dt;
                
                
                
                
                
            end
            
            
            
            
            
            allOTe(aii,:,:)=[Re]; %; Rpv; Rvip; Rsom];
            allOTpv(aii,:,:)=[Rpv];
            allOTvip(aii,:,:)=[Rvip];
            allOTsom(aii,:,:)=[Rsom];
        end
        
        
        %peak:
        OTeP(ti,:,:)=squeeze(max(allOTe(:,:,endT-stT:endT),[],3));
        OTpvP(ti,:,:)=squeeze(max(allOTpv(:,:,endT-stT:endT),[],3));
        OTvipP(ti,:,:)=squeeze(max(allOTvip(:,:,endT-stT:endT),[],3));
        OTsomP(ti,:,:)=squeeze(max(allOTsom(:,:,endT-stT:endT),[],3));
        %mean transient:
        OTeM(ti,:,:)=squeeze(mean(allOTe(:,:,stT:stT+100),3));
        OTpvM(ti,:,:)=squeeze(mean(allOTpv(:,:,stT:stT+100),3));
        OTvipM(ti,:,:)=squeeze(mean(allOTvip(:,:,stT:stT+100),3));
        OTsomM(ti,:,:)=squeeze(mean(allOTsom(:,:,stT:stT+100),3));
        %steady state:
        OTeT(ti,:,:)=squeeze(mean(allOTe(:,:,endT-10:endT),3));
        OTpvT(ti,:,:)=squeeze(mean(allOTpv(:,:,endT-10:endT),3));
        OTvipT(ti,:,:)=squeeze(mean(allOTvip(:,:,endT-10:endT),3));
        OTsomT(ti,:,:)=squeeze(mean(allOTsom(:,:,endT-10:endT),3));
        %baseline:
        OTeTs(ti,:,:)=squeeze(mean(allOTe(:,:,200:499),3));
        OTpvTs(ti,:,:)=squeeze(mean(allOTpv(:,:,200:499),3));
        OTvipTs(ti,:,:)=squeeze(mean(allOTvip(:,:,200:499),3));
        OTsomTs(ti,:,:)=squeeze(mean(allOTsom(:,:,200:499),3));
        
    end
    
    
    
    mE=squeeze(mean(OTeT,1)); mP=squeeze(mean(OTpvT,1)); mV=squeeze(mean(OTvipT,1)); mS=squeeze(mean(OTsomT,1));
    mEs=squeeze(mean(OTeTs,1)); mPs=squeeze(mean(OTpvTs,1)); mVs=squeeze(mean(OTvipTs,1)); mSs=squeeze(mean(OTsomTs,1));
    
    
    aveEv(ni,:)=[mean(mean(OTeT(:,2,:),3),1), mean(mean(OTpvT(:,2,:),3),1), mean(mean(OTvipT(:,2,:),3),1), mean(mean(OTsomT(:,2,:),3),1)];
    
    aveSp(ni,:)=[mean(mean(OTeTs(:,2,:),3),1), mean(mean(OTpvTs(:,2,:),3),1), mean(mean(OTvipTs(:,2,:),3),1), mean(mean(OTsomTs(:,2,:),3),1)];
    
    Be=[mE(1,:)-mE(2,:)  ];
    Bes=[mEs(1,:)-mEs(2,:)  ];
    
    
    Bp=[mP(1,:)-mP(2,:)  ];
    Bps=[mPs(1,:)-mPs(2,:)  ];
    
    
    Bv=[mV(1,:)-mV(2,:)  ];
    Bvs=[mVs(1,:)-mVs(2,:)  ];
    
    
    Bs=[mS(1,:)-mS(2,:)  ];
    Bss=[mSs(1,:)-mSs(2,:)  ];
    
    avedel(ni,:)=[sum(Be)/length(Be), sum(Bp)/length(Bp), sum(Bv)/length(Bv), sum(Bs)/length(Bs)];
    avedels(ni,:)=[sum(Bes)/length(Bes), sum(Bps)/length(Bps), sum(Bvs)/length(Bvs), sum(Bss)/length(Bss)];
    
    
    aveperc(ni,:)=[sum(Be<0)/length(Be), sum(Bp>0)/length(Bp), sum(Bv<0)/length(Bv), sum(Bs>0)/length(Bs)];
    avepercs(ni,:)=[sum(Bes>0)/length(Bes), sum(Bps>0)/length(Bps), sum(Bvs>0)/length(Bvs), sum(Bss>0)/length(Bss)];
    
    
    
end








f3=figure
subplot(2,3,1)
hold on
bar(mean(aveSp,1))
errorbar(mean(aveSp,1),std(aveSp,[],1),'.k')
SpontRate(1,:)=mean(aveSp,1); SpontRate(2,:)=std(aveSp,[],1);
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
title('Spontaneous Firing')
subplot(2,3,4)
hold on
bar(mean(aveEv,1))
EvRatesSS(1,:)=mean(aveEv,1); EvRatesSS(2,:)=std(aveEv,[],1);
errorbar(mean(aveEv,1),std(aveEv,[],1),'.k')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
title('Evoked Firing')
subplot(2,3,2)
hold on
bar(mean(avedels,1))
spDifs(1,:)=mean(avedels,1); spDifs(2,:)=std(avedels,[],1);
errorbar(mean(avedels,1),std(avedels,[],1),'.k')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
title('Ach Rate - No Ach')

subplot(2,3,3)
hold on
bar(mean(avepercs,1))
spPerc(1,:)=mean(avepercs,1) ; spPerc(2,:)= std(avepercs,[],1);
errorbar(mean(avepercs,1),std(avepercs,[],1),'.k')
%scatter([1 2 3 4],[.17,.23, .28, .1]+.5, 'r')
plot([1 4], [.5 .5])
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})

title('% cells increased')

subplot(2,3,5)
hold on
bar(mean(avedel,1))
DifsSS(1,:)=mean(avedel,1); DifsSS(2,:)=std(avedel,[],1);
errorbar(mean(avedel,1),std(avedel,[],1),'.k')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
title('Ach Rate - No Ach')

subplot(2,3,6)
hold on
bar(mean(aveperc,1))
PercSS(1,:)=mean(aveperc,1); PercSS(2,:)=std(aveperc,[],1);
errorbar(mean(aveperc,1),std(aveperc,[],1),'.k')
scatter([1 2 3 4],[.17,.11, .28, .1]+.5, 'r')
plot([1 4], [.5 .5])
ylim([0 1])
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})

title('% cells changed in the right direction')









