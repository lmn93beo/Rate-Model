rng('shuffle');

NumNeurons = 400;
time = linspace(0,0.8,10000);
for sim = 1:10
    
    for neuron = 1:NumNeurons
        
        mu =  0.0 + 0.75*rand(1,1);
        sigma =  0.1 + 0.01*rand(1,1);
        
        MU(neuron) = mu; SIG(neuron) = sigma;
        
        pd = makedist('Normal',mu,sigma);
        ModelPopA(neuron,:) = pdf(pd, time) + normrnd(0.02, 0.05,1,length(time));
        [amp(neuron),m(neuron)]  = max(ModelPopA(neuron,:));
        
        ModelPopB(neuron,:) = 1*normrnd(mu, sigma ,1,length(time));
        
    end;
    [v,ind] = sort(m);
    figure(1);
    subplot(1,3,1); imagesc(time,1:NumNeurons, ModelPopA(ind,:) );
    subplot(1,3,2); imagesc(time,1:NumNeurons, ModelPopB(:,:) );
    
    %%
    
    Net = [ModelPopA;ModelPopB];
    Net = Net( randperm(size(Net,1)), :);
    
    figure(1);
    subplot(1,3,3); imagesc(time, 1:size(Net,1), Net );
    
    [coeff, score, latent, tsq, expl, mue] = pca(Net);
    ExplVar = sum( expl(1:3) );
    figure(2);
    plot3( coeff(:,1), coeff(:,2), coeff(:,3) );
    hold on;
    
    plot3( coeff(1,1), coeff(1,2), coeff(1,3) ,'o', 'markersize', 12,'markerfacecolor','r');
    plot3( coeff(end,1), coeff(end,2), coeff(end,3) ,'o', 'markersize', 12,'markerfacecolor','g');
    
    disp( ['Sim #: ' num2str(sim) '  Expl.Var: ' num2str(ExplVar)] );
    
    
end;
%%
proj = score(:,1:10)*coeff(:,1:10)';
proj = proj(1:10,:);

figure(3);
subplot(1,3,1); plot(cumsum(expl(1:10)),'-ro'); hold on;
subplot(1,3,2); imagesc(time, 1:10, proj );
subplot(1,3,3); histogram(MU, 40);


