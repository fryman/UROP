function LRatio = af_rateClusters(s)
% Reads cluster feature data from .MAT file <fn1>, and returns a list of
% quality ratings with one element for each cluster in the file.  Cluster 0
% is not included in <LRatio>, but it must be supplied in <fn1> in order to
% do the computation.  The quality computation is the L(ratio) described in
% the CLUSTER QUALITY section of Schmitzer-Torbert N & Redish AD (2004) J
% Neurophysiol 91:2259-2272. For explanation of Mahalanobis distance, see
% Fall 2005 web notes for Princeton Computer Science 436 at
% http://www.cs.princeton.edu/courses/archive/fall05/cos436/Duda/PR_Mahal/M_metric.htm. 
%INPUT
% The file specified by <fn1> is a .MAT file containing cluster ID numbers
% in column 2 and spike feature values in subsequent columns.  Warning:
% computation time goes as the square of the number of features.  Column 1
% is ignored.  The file must contain "unassigned" or "noise" samples as
% cluster #0 in order to get a legitimate assessment of quality.  The
% cluster numbers in <fn1> do not need to be consecutive; any missing
% cluster numbers are simply skipped.

%$Rev: 57 $
%$Date: 2010-06-19 14:05:00 -0400 (Sat, 19 Jun 2010) $
%$Author: dgibson $

% Load the data into arrays s:
% S = load(fn1);
% fields = fieldnames(S);
% s = S.(fields{1});
% clear S;

clusters = unique(s(:,2))';
if ~ismember(0, clusters)
    error('dg_rateClusters:no0', ...
        'There is no cluster #0 in %s', s );
end
Nfeatures = size(s, 2) - 2;
if Nfeatures < 1
    error('dg_rateClusters:nofeat', ...
        'There are no feature data in %s', s );
end
disp(sprintf('Rating clusters based on %d features', Nfeatures));
for clustnum = clusters(2:end)
    clustselect = (s(:,2) == clustnum);
    % I make the leap of faith that the incredibly compact formulas in the
    % Matlab 'mahal' function really do compute the Mahalanobis distances
    % of Y from mean(X) using the the covariance matrix calculated from X.
    L(clustnum) = sum(1 - chi2cdf(...
        mahal(s(clustselect, 3:end), s(~clustselect, 3:end)), ...
        Nfeatures ));
    LRatio(clustnum) = L(clustnum) / sum(clustselect);
end