% function rateModelAccum


%% Toy rate model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rng('shuffle')

%% Connection Weights >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Wep =  to e from
W11 = 0.987;    W12 = 0;        W13 = 0.512;   W14 = 0;       W15 = -0.68;     W16 = 0; % inputs to epop1 (pfc 1)
W21 = 0;        W22 = 0.987;    W23 = 0.512;   W24 = -0.68;   W25 = 0;         W26 = 0;   % inputs to epop2 (pfc 2)
W31 = 1.385;    W32 = 1.385;    W33 = 0;        W34 = 0;       W35 = 0;         W36 = 0;    % inputs to epop3 (thal)
W41 = 0.685;    W42 = 0;        W43 = 0.53;     W44 = -0.12;   W45 = 0;         W46 = 0;    % inputs to ipop3 (inh 1)
W51 = 0;        W52 = 0.685;    W53 = 0.53;     W54 = 0;       W55 = -0.12;     W56 = 0;  % inputs to ipop3 (inh 1)
W61 = 1.385;    W62 = 0 ;       W63 = 0.00;     W64 = 0;       W65 = 0;         W66 = 0; % inputs to ipop4 (trn)

f = @(I, k, n) k.*( ( I.*(I>0) ).^n );

x = 1;

%% define time constants...................................................
tau_1 = 25;
tau_2 = 25;
tau_3 = 20;
tau_4 = 10;
tau_5 = 10;
tau_6 = 10;

DelT  = 0.05; %integration time step in ms.
time  = -200:DelT:900;


%% define simulation parameters
NumTrials = 60;

isPulseStim      = 0;
isVariableWeight = 1;
isStimNoisy      = 1;
isBGNoisy        = 1;
diagnosticMode   = 0;
anaCurrents      = 0;


plotFlag         = 0;

Istim1_amp =  10;
Istim2_amp =  10;


ind = find( time >= 400 & time <= 600 );
gain = ones( 1, length(time) );
y  = exp(-0.08*time(ind));
gain(ind) = y./max(y);

%%

for stim = 1:10
    vispure = [1,0,0,0,0,0,1,0,0,0];
    audpure = [1,0,1,1,1,1,1,0,1,1];
    
    vismix = [1,0,0,1,0,1,0,1,0,0];
    audmix = [1,0,1,0,1,0,1,1,0,1];
    
    vismix = vismix(randperm(10));
    audmix = audmix(randperm(10));
    
    
    for trial = 1:NumTrials
        % Define stimulus >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        
        
        Istim1 = zeros(1,length(time));
        Istim2 = zeros(1,length(time));
        
        durPulse = 50; durISI = 25;
        Npulse = 9;
        
        tvec = [0:(durPulse+durISI):Npulse*(durPulse+durISI)];
        
        for i = 1:length(tvec)
            ind1 = find( time >= tvec(i) & time <= durPulse+tvec(i) );
            Istim1(ind1) = vismix(i).*Istim1_amp.*ones(1,length(ind1));
            Istim2(ind1) = audmix(i).*Istim2_amp.*ones(1,length(ind1));
        end;
        
        Istim1 = smoothts(Istim1,'e',500);
        Istim2 = smoothts(Istim2,'e',500);
        
        if mod(trial,2) == 1
            W1s = 0; W2s = 0; W3s = 1;
            W4s = 1; W5s = 0;
        else
            W1s = 0; W2s = 0; W3s = 1;
            W4s = 0; W5s = 1;
        end;
        
        
        
        clear W;
        
        W = [ W11, W12, W13, W14, W15, W16;...
            W21, W22, W23, W24, W25, W26;...
            W31, W32, W33, W34, W35, W36; ...
            W41, W42, W43, W44, W45, W46; ...
            W51, W52, W53, W54, W55, W56; ...
            W61, W62, W63, W64, W65, W66];
        
        if isVariableWeight
            
            W = W.*( 0.9+0.2*rand(size(W)) );
            Weight{trial} = W;
            
        end;
        
        % initialise model.
        r_1  = 1.*ones(1,length(time));
        r_2  = 1.*ones(1,length(time));
        r_3  = 1.*ones(1,length(time));
        r_4 =  8.*ones(1,length(time));
        r_5 =  8.*ones(1,length(time));
        r_6 =  8.*ones(1,length(time));
        
        % background currents
        I1_b = poissrnd(10.8,1,length(time));
        I2_b = poissrnd(10.8,1,length(time));
        I3_b = poissrnd(15,1,length(time));
        I4_b = poissrnd(20,1,length(time));
        I5_b = poissrnd(20,1,length(time));
        I6_b = poissrnd(20,1,length(time));
        
        for t = 1:length(time)-1
            
            
            
            I1 = W(1,1).*r_1(t) + W(1,2).*r_2(t) + W(1,3).*r_3(t) + W(1,4)*r_4(t) + W(1,5)*r_5(t)+ W(1,6)*r_6(t)+  W1s*Istim1(t) + I1_b(t);
            I2 = W(2,1).*r_1(t) + W(2,2).*r_2(t) + W(2,3).*r_3(t) + W(2,4)*r_4(t) + W(2,5)*r_5(t)+ W(2,6)*r_6(t)+  W2s*Istim2(t) + I2_b(t);
            I3 = W(3,1).*r_1(t) + W(3,2).*r_2(t) + W(3,3).*r_3(t) + W(3,4)*r_4(t) + W(3,5)*r_5(t)+ W(3,6)*r_6(t)+  I3_b(t) + W3s*Istim2(t) + W3s*Istim1(t);
            I4 = W(4,1).*r_1(t) + W(4,2).*r_2(t) + W(4,3).*r_3(t) + W(4,4)*r_4(t) + W(4,5)*r_5(t)+ W(4,6)*r_6(t)+  I4_b(t) + W4s*Istim1(t);
            I5 = W(5,1).*r_1(t) + W(5,2).*r_2(t) + W(5,3).*r_3(t) + W(5,4)*r_4(t) + W(5,5)*r_5(t)+ W(5,6)*r_6(t)+  I5_b(t) + W5s*Istim2(t);
            I6 = W(6,1).*r_1(t) + W(6,2).*r_2(t) + W(6,3).*r_3(t) + W(6,4)*r_4(t) + W(6,5)*r_5(t)+ W(6,6)*r_6(t)+  I6_b(t);
            
            r_1(t+1) = r_1(t) + (DelT/tau_1).*(-r_1(t) + f(I1,0.01,2.2) );
            r_2(t+1) = r_2(t) + (DelT/tau_2).*(-r_2(t) + f(I2,0.01,2.2) );
            r_3(t+1) = r_3(t) + (DelT/tau_3).*(-r_3(t) + f(I3,0.01,2.2) );
            r_4(t+1) = r_4(t) + (DelT/tau_4).*(-r_4(t) + f(I4,0.01,2.2) );
            r_5(t+1) = r_5(t) + (DelT/tau_5).*(-r_5(t) + f(I5,0.01,2.2) );
            r_6(t+1) = r_6(t) + (DelT/tau_6).*(-r_6(t) + f(I6,0.01,2.2) );
        end;
        
        Rate{1,trial} = ( r_1 - min(r_1) )./( max(r_1)-min(r_1));
        Rate{2,trial} = ( r_2 - min(r_2) )./( max(r_2)-min(r_2));
        Rate{3,trial} = ( r_3 - min(r_3) )./( max(r_3)-min(r_3));
        Rate{4,trial} = ( r_4 - min(r_4) )./( max(r_4)-min(r_4));
        Rate{5,trial} = ( r_5 - min(r_5) )./( max(r_5)-min(r_5));
        Rate{6,trial} = ( r_6 - min(r_6) )./( max(r_6)-min(r_6));
        
        %%
        if mod(trial, 2) == 1
            figure(1); set(gcf,'color','w');
            
            subplot(1,3,1)
            plot( time, r_1, 'r'); hold on;
            subplot(1,3,2);
            plot( time, r_2, 'b'); hold on;
            subplot(1,3,3);
            plot( time, r_3, 'k'); hold on;
        else
            
            figure(2); set(gcf,'color','w');
            
            subplot(1,3,1)
            plot( time, r_1, 'r'); hold on;
            subplot(1,3,2);
            plot( time, r_2, 'b'); hold on;
            subplot(1,3,3);
            plot( time, r_3, 'k'); hold on;
        end;
        
    end;
    
    
    %%
    rr1 = cell2mat( Rate(1,:)' );
    rr2 = cell2mat( Rate(2,:)' );
    rr3 = cell2mat( Rate(3,:)' );
    
    figure(3);
    % plot( time, mean(rr2(1:2:end, :)),'color','r', 'linewidth',3); hold on;
    plot( time, mean(rr1(1:2:end, :)),'color','b', 'linewidth',3); hold on;
    plot( time, mean(rr3(1:2:end, :)),'color','k', 'linewidth',3); hold on;
    
    figure(4);
    % plot( time, mean(rr1(2:2:end, :)),'color','r', 'linewidth',3); hold on;
    plot( time, mean(rr2(2:2:end, :)),'color','b', 'linewidth',3); hold on;
    plot( time, mean(rr3(2:2:end, :)),'color','k', 'linewidth',3); hold on;
end
