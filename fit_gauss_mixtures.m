function [optK,bestpp,optMu,optCovar,dl,countf,optIndicators,trackK,bestNormIndicators] = fit_gauss_mixtures(y, w, kmin, kmax, regularize, th, weight_threshold, rand_indices)
%
%
% Usage syntax:
% [optK,bestpp,optMu,optCovar,dl,countf] = mixtures3(y,kmin,kmax,regularize,th,covoption)
%
% Inputs:
%
% "y" is the data; for n observations in d dimensions, y must have
% d lines and n columns.
%
% "kmax" is the initial (maximum) number of mixture components
% "kmin" is the minimum number of components to be tested
%
% "regularize" is a regularizing factor for covariance matrices; in very small samples,
% it may be necessary to add some small quantity to the diagonal of the covariances
%
% "th" is a stopping threshold
%
% Outputs:
%
% "optK" is the selected number of components
% "bestpp" is the obtained vector of mixture probabilities
% "optMu" contains the estimates of the means of the components
%          it has optK columns by d lines
% "optCovar" contains the estimates of the covariances of the components
%           it is a three-dimensional (three indexes) array
%           such that optCovar(:,:,1) is the d by d covariance of the first
%           component and optCovar(:,:,optK) is the covariance of the last
%           component
%
% "dl" contains the successive values of the cost function
% "countf" is the total number of iterations performed

verb=0; % verbose mode; change to zero for silent mode
[dimens,npoints] = size(y);
if w == 1
    w = ones(1,size(y,2));
end


npars = (dimens + dimens*(dimens+1)/2);

nparsover2 = npars / 2;

weights = w*npoints/sum(w);

% kmax is the initial number of mixture components
k = kmax;

% indic variables are indicator functions and will
% contain the assignments of each data point to
% the mixture components, as result of the E-step
indic = zeros(k,npoints);

% Initialization: we will initialize the means of the k components
% with k randomly chosen data points. Randperm(n) is a MATLAB function
% that generates random permutations of the integers from 1 to n.

% indices = find(w>weight_threshold);
% randindex = randperm(size(indices,2));
% randindex = indices(randindex(1:k));

% randindex = randperm(npoints);
% randindex = randindex(1:k);

%estMu is the current best estimate of the means
% estMu = y(:,randindex);
% rand_indices = ceil(npoints*rand_indices);
estMu = y(:,rand_indices);

% the initial estimates of the mixing probabilities are set to 1/k
estpp = (1/k)*ones(1,k);

% here we compute the global covariance of the data
% globCovar = cov(y');  %DBUG--WEIGHTS

% the covariances are initialized to 1/10 the global covariance.
% this is a completely aribitrary choce and can be changed. Original
% code set covariance to an isotropic gaussian with diagonal covariance
% for most cases the newer version converges faster.
estCovar(:,:,1:k) = repmat(0.1*cov(y'), [1,1,k]);

% having the initial means, covariances, and probabilities, we can
% initialize the indicator functions following the standard EM equation
% Notice that these are unnormalized values
%keyboard
for i=1:k
    semi_indic(i,:) = multinorm(y,estMu(:,i),estCovar(:,:,i));
end

indic(1:k,:) = semi_indic(1:k,:).*repmat(estpp(1:k)', [1,length(semi_indic)]);


% we can use the indic variables (unnormalized) to compute the
% loglikelihood and store it for later plotting its evolution
% we also compute and store the number of components
countf = 1;

loglike(countf) = sum(weights.*log(sum(realmin+indic)));
dlength = -loglike(countf) + (nparsover2*sum(log(estpp))) + ...
    (nparsover2 + 0.5)*k*log(sum(weights));
dl(countf) = dlength;

% the transitions vectors will store the iteration
% number at which components are killed.
% transitions1 stores the iterations at which components are
% killed by the M-step, while transitions2 stores the iterations
% at which we force components to zero.

% minimum description length seen so far, and corresponding
% parameter estimates
minDL = dl(countf);
bestpp = estpp;
optMu = estMu;
optCovar = estCovar;
optK = k;
optIndicators = indic;
bestNormIndicators = indic./(realmin+repmat(sum(indic,1), [k,1])).*repmat(weights,[size(indic,1),1]);

transitions1 = [];
transitions2 = [];

outer = 0;

disp('Clustering...')
while outer == 0
    inner = 0;
    while inner == 0
        % this inner loop is the component-wise EM algorithm with the
        % modified M-step that can kill components
        
        if verb==1
            % in verbose mode, we keep displaying the minimum of the
            % mixing probability estimates to see how close we are
            % to killing one component
            fprintf('k = %2d,  minestpp = %0.5g', k, min(estpp));
        end
        
        comp = 1;
        % Since k may change during the process, we can not use a for loop
        while comp <= k
            % we start with the M step
            % first, we compute a normalized indicator function
            clear indic
            indic(1:k,:) = semi_indic(1:k,:).*repmat(estpp(1:k)', ...
                [1,length(semi_indic)]);
            normindic = indic./(realmin+repmat(sum(indic,1), [k,1])).*repmat(weights,[size(indic,1),1]);
            
            % now we perform the standard M-step for mean and covariance
            normalize = 1/sum(normindic(comp,:));
            weightedY = repmat(normindic(comp,:),[dimens,1]).*y;
            estMu(:,comp)= normalize*sum(weightedY,2);
            estCovar(:,:,comp) = normalize*(weightedY*y') - estMu(:,comp)*estMu(:,comp)' ...
                + regularize*eye(dimens);
            
            
            % this is the special part of the M step that is able to
            % kill components
            compute_probs();
            
            % this is an auxiliary variable that will be used the
            % signal the killing of the current component being updated
            killed = 0;
            
            % we now have to do some book-keeping if the current component was killed
            % that is, we have to rearrange the vectors and matrices that store the
            % parameter estimates
            if estpp(comp)==0
                killed = 1;
                % we also register that at the current iteration a component was killed
                transitions1 = [transitions1 countf];
                
                estMu = [estMu(:,1:comp-1) estMu(:,comp+1:k)];
                newcov = zeros(dimens,dimens,k-1);
                for kk=1:comp-1
                    newcov(:,:,kk) = estCovar(:,:,kk);
                end
                for kk=comp+1:k
                    newcov(:,:,kk-1) = estCovar(:,:,kk);
                end
                estCovar = newcov;
                estpp = [estpp(1:comp-1) estpp(comp+1:k)];
                semi_indic = semi_indic([1:comp-1,comp+1:k],:);
                
                % since we've just killed a component, k must decrease
                k = k-1;
            end
            
            if killed==0
                % if the component was not killed, we update the corresponding
                % indicator variables...
                semi_indic(comp,:) = multinorm(y,estMu(:,comp),estCovar(:,:,comp));
                % ...and go on to the next component
                comp = comp + 1;
            end
            % if killed==1, it means the estCovarin the position "comp", we now
            % have a component that was not yet visited in this sweep, and
            % so all we have to do is go back to the M setp without
            % increasing "comp"
            
        end
        countf = countf + 1;
        compute_dl();
%         display(size(dl));
        % compute the change in loglikelihood to check if we should stop
        deltlike = loglike(countf) - loglike(countf-1);
        if (verb~=0)
            fprintf('deltaloglike/th = %0.7g', abs(deltlike/loglike(countf-1))/th);
        end
        if (abs(deltlike/loglike(countf-1)) < th)
            % if the relative change in loglikelihood is below the threshold, we stop CEM2
            inner = 1;
        end
        
    end
       
    % now check if the latest description length is the best;
    % if it is, we store its value and the corresponding estimates
    if dl(countf) < minDL
        bestpp = estpp;
        optMu = estMu;
        optCovar = estCovar;
        optK = k;
        minDL = dl(countf);
        optIndicators = indic;
        bestNormIndicators = indic./(realmin+repmat(sum(indic,1), [k,1])).*repmat(weights,[size(indic,1),1]);
    end
    
    % at this point, we may try smaller mixtures by killing the
    % component with the smallest mixing probability and then restarting CEM2,
    % as long as k is not yet at kmin
    if k > kmin
        [~, indminp] = min(estpp);
        % what follows is the book-keeping associated with removing one component
        if indminp==1
            estMu = estMu(:,2:k);
            estCovar = estCovar(:,:,2:k);
            estpp = estpp(2:k);
        elseif indminp==k
            estMu = estMu(:,1:k-1);
            estCovar = estCovar(:,:,1:k-1);
            estpp = estpp(1:k-1);
        else
            estMu = [estMu(:,1:indminp-1) estMu(:,indminp+1:k)];
            newcov = zeros(dimens,dimens,k-1);
            for kk=1:indminp-1
                newcov(:,:,kk) = estCovar(:,:,kk);
            end
            for kk=indminp+1:k
                newcov(:,:,kk-1) = estCovar(:,:,kk);
            end
            estCovar = newcov;
            estpp = [estpp(1:indminp-1) estpp(indminp+1:k)];
        end
        k = k-1;
        
        % we renormalize the mixing probabilities after killing the component
        estpp = estpp/sum(estpp);
        
        % and register the fact that we have forced one component to zero
        transitions2 = [transitions2 countf];
        
        % ...and compute the loglikelihhod function and the description length
        countf = countf + 1;
        compute_dl();
    else
        outer = 1;
    end
end

    function compute_dl()
        
        %       clear semi_indic; clear indic;
        indic = calc_indicators(y,estMu,estCovar,k,estpp);
        %         for i=1:k
        %             semi_indic(i,:) = multinorm(y,estMu(:,i),estCovar(:,:,i));
        %         end
        %
        %         indic(1:k,:) = semi_indic(1:k,:).*repmat(estpp(1:k)', [1,length(semi_indic)]);
        
        if k~=1
            % if the number of surviving components is not just one, we compute
            % the loglikelihood from the unnormalized assignment variables
            loglike(countf) = sum(weights.*log(realmin+sum(indic)./weights)); %.*repmat(W,[size(indic,1),1]))));
        else
            % if it is just one component, it is even simpler
            loglike(countf) = sum(weights.*log(realmin+indic./weights)); %.*repmat(W,[size(indic,1),1])));
        end
        
        % compute and store the description length and the current number of components
        dl(countf) = -loglike(countf) + (nparsover2*sum(log(estpp))) + ...
            (nparsover2 + 0.5)*k*log(sum(weights));
        trackK(countf) = k;
    end


    function compute_probs()
        estpp(comp) = max(sum(normindic(comp,:))-nparsover2,0)/sum(weights);
        estpp = estpp/sum(estpp);
    end

    function [indicators] = calc_indicators(data, means, covars, comps, mixProbs)
        for idx=1:comps
            semi_indic(idx,:) = multinorm(data,means(:,idx),covars(:,:,idx));
        end
        
        indicators(1:comps,:) = semi_indic(1:comps,:).*repmat(mixProbs(1:comps)', [1,length(semi_indic)]);
    end

end



