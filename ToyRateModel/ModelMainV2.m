
%% Toy rate model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rng('shuffle')

%% Connection Weights >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Wep =  to e from
W11 = 0.987;    W12 = 0;        W13 = 0.412;   W14 = 0;       W15 = -0.65;     W16 = 0; % inputs to epop1 (pfc 1)
W21 = 0;        W22 = 0.987;    W23 = 0.412;   W24 = -0.65;   W25 = 0;         W26 = 0;   % inputs to epop2 (pfc 2)
W31 = 1.385;    W32 = 1.385;    W33 = 0;        W34 = 0;       W35 = 0;         W36 = -1.995;    % inputs to epop3 (thal)
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

W1s = 1; W2s = 1; W3s = 1; 
W4s = 0; W5s = 0;

ind = find( time >= 400 & time <= 600 );
gain = ones( 1, length(time) );
y  = exp(-0.08*time(ind));
gain(ind) = y./max(y);

%%

for trial = 1:NumTrials
    % Define stimulus >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    trialType = randi(4);
    disp(['Trial Type: ' num2str(trialType)]);
    
    if trialType ==  1
        % single R1
        
        Istim1 = zeros(1,length(time));
        Istim2 = zeros(1,length(time));
        
        tstart1 = normrnd(50, 45);
        tstart2 = normrnd(550, 45);
        
        ind1 = find( time >= tstart1 & time <= 100+tstart1 );
        ind2 = find( time >= tstart2 & time <= 100+tstart2 );
        
        Y = genrateGamCurrent(100);
        
        Istim1( ind1  )  =  Istim1_amp.*Y(1:end-1);
        
    elseif trialType == 2
        % single R2
        Istim1 = zeros(1,length(time));
        Istim2 = zeros(1,length(time));
        
        tstart1 = normrnd(50, 45);
        tstart2 = normrnd(550, 45);
        
        ind1 = find( time >= tstart1 & time <= 100+tstart1 );
        ind2 = find( time >= tstart2 & time <= 100+tstart2 );
        
        Y = genrateGamCurrent(100);
        
        Istim2( ind1  )  =  Istim2_amp.*Y(1:end-1);
        
    elseif trialType == 3
        % repeated R1
        
        Istim1 = zeros(1,length(time));
        Istim2 = zeros(1,length(time));
        
        tstart1 = normrnd(50, 45);
        tstart2 = normrnd(550, 45);
        
        ind1 = find( time >= tstart1 & time <= 100+tstart1 );
        ind2 = find( time >= tstart2 & time <= 100+tstart2 );
        
        Y = genrateGamCurrent(100);
        
        Istim1( ind1  )  =  Istim1_amp.*Y(1:end-1);
        Istim1( ind2  )  =  Istim1_amp.*Y(1:end-1);
        
    elseif trialType == 4
        % conflict
        Istim1 = zeros(1,length(time));
        Istim2 = zeros(1,length(time));
        
        tstart1 = normrnd(50, 45);
        tstart2 = normrnd(550, 45);
        
        ind1 = find( time >= tstart1 & time <= 100+tstart1 );
        ind2 = find( time >= tstart2 & time <= 100+tstart2 );
        
        Y = genrateGamCurrent(100);
        
        Istim1( ind1  )  =  Istim1_amp.*Y(1:end-1);
        Istim2( ind2  )  =  Istim2_amp.*Y(1:end-1);
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
        
        if trialType == 3 || trialType == 4;
            
            W(3,6)= gain(t).*W36;
        end;
        
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
    
    Rate{1,trial} = r_1;
    Rate{2,trial} = r_2;
    Rate{3,trial} = r_3;
    Rate{4,trial} = r_4;
    Rate{5,trial} = r_5;
    Rate{6,trial} = r_6;
    
    %%
    figure(1); set(gcf,'color','w');
    if trialType == 1
        subplot(1,4,1);
        plot( time, r_3, 'k'); hold on;
         ylim([0,30]);
    elseif trialType == 2
        subplot(1,4,2);
        plot( time, r_3, 'k'); hold on;
         ylim([0,30]);
    elseif trialType == 3
        subplot(1,4,3);
        plot( time, r_3, 'k'); hold on;
         ylim([0,30]);
    elseif trialType == 4
        subplot(1,4,4);
        plot( time, r_3, 'k'); hold on;
         ylim([0,30]);
    end;
    
    figure(2); set(gcf,'color','w');
    if trialType == 1
        subplot(1,4,1);
        plot( time, r_1, 'r'); hold on;
        plot( time, r_2, 'b'); hold on;
%         ylim([0,50]);
        
    elseif trialType == 2
        subplot(1,4,2);
        plot( time, r_1, 'r'); hold on;
        plot( time, r_2, 'b'); hold on;
%         ylim([0,50]);
    elseif trialType == 3
        subplot(1,4,3);
        plot( time, r_1, 'r'); hold on;
        plot( time, r_2, 'b'); hold on;
%         ylim([0,50]);
    elseif trialType == 4
        subplot(1,4,4);
        plot( time, r_1, 'r'); hold on;
        plot( time, r_2, 'b'); hold on;
%         ylim([0,50]);
    end;
    
    
    figure(3); set(gcf,'color','w');
    if trialType == 1
        subplot(1,4,1);
        plot( time, r_6, 'r'); hold on;
%         ylim([0,50]);
        
    elseif trialType == 2
        subplot(1,4,2);
        plot( time, r_6, 'r'); hold on;
%         ylim([0,50]);
    elseif trialType == 3
        subplot(1,4,3);
        plot( time, r_6, 'r'); hold on;
%         ylim([0,50]);
    elseif trialType == 4
        subplot(1,4,4);
        plot( time, r_6, 'r'); hold on;
%         ylim([0,50]);
    end;
    
end;