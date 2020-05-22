function varargout = uniform_kmeans(X, k, varargin)
% couturier: this a small modification of the Matlab kmeans function. It forces clusters of uniform size. The modification is in the uniform_min function at the end.
%KMEANS K-means clustering.
%   IDX = KMEANS(X, K) partitions the points in the N-by-P data matrix X
%   into K clusters.  This partition minimizes the sum, over all clusters, of
%   the within-cluster sums of point-to-cluster-centroid distances.  Rows of X
%   correspond to points, columns correspond to variables.  Note: when X is a
%   vector, KMEANS treats it as an N-by-1 data matrix, regardless of its
%   orientation.  KMEANS returns an N-by-1 vector IDX containing the cluster
%   indices of each point.  By default, KMEANS uses squared Euclidean
%   distances.
%
%   KMEANS treats NaNs as missing data, and ignores any rows of X that
%   contain NaNs.
%
%   [IDX, C] = KMEANS(X, K) returns the K cluster centroid locations in
%   the K-by-P matrix C.
%
%   [IDX, C, SUMD] = KMEANS(X, K) returns the within-cluster sums of
%   point-to-centroid distances in the K-by-1 vector sumD.
%
%   [IDX, C, SUMD, D] = KMEANS(X, K) returns distances from each point
%   to every centroid in the N-by-K matrix D.
%
%   [ ... ] = KMEANS(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by KMEANS.  Parameters are:
%
%   'Distance' - Distance measure, in P-dimensional space, that KMEANS
%      should minimize with respect to.  Choices are:
%          'sqeuclidean'  - Squared Euclidean distance (the default)
%          'cityblock'    - Sum of absolute differences, a.k.a. L1 distance
%          'cosine'       - One minus the cosine of the included angle
%                           between points (treated as vectors)
%          'correlation'  - One minus the sample correlation between points
%                           (treated as sequences of values)
%          'hamming'      - Percentage of bits that differ (only suitable
%                           for binary data)
%
%   'Start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%          'plus'    - The Default. Select K observations from X according
%                      to the k-means++ algorithm: the first cluster center
%                      is chosen uniformly at random from X, after which
%                      each subsequent cluster center is chosen randomly
%                      from the remaining data points with probability
%                      proportional to its distance from the point's
%                      closest existing cluster center.
%          'sample'  - Select K observations from X at random.
%          'uniform' - Select K points uniformly at random from the range
%                      of X.  Not valid for Hamming distance.
%          'cluster' - Perform preliminary clustering phase on random 10%
%                      subsample of X.  This preliminary phase is itself
%                      initialized using 'sample'.
%           matrix   - A K-by-P matrix of starting locations.  In this case,
%                      you can pass in [] for K, and KMEANS infers K from
%                      the first dimension of the matrix.  You can also
%                      supply a 3D array, implying a value for 'Replicates'
%                      from the array's third dimension.
%
%   'Replicates' - Number of times to repeat the clustering, each with a
%      new set of initial centroids.  A positive integer, default is 1.
%
%   'EmptyAction' - Action to take if a cluster loses all of its member
%      observations.  Choices are:
%          'singleton' - Create a new cluster consisting of the one
%                        observation furthest from its centroid (the
%                        default).
%          'error'     - Treat an empty cluster as an error.
%          'drop'      - Remove any clusters that become empty, and set
%                        the corresponding values in C and D to NaN.
%          
%
%   'Options' - Options for the iterative algorithm used to minimize the
%       fitting criterion, as created by STATSET.  Choices of STATSET
%       parameters are:
%
%          'Display'  - Level of display output.  Choices are 'off', (the
%                       default), 'iter', and 'final'.
%          'MaxIter'  - Maximum number of iterations allowed.  Default is 100.
%
%      'UseParallel'  - If true and if a parpool of the Parallel Computing
%                       Toolbox is open, compute in parallel. If the
%                       Parallel Computing Toolbox is not installed, or a
%                       parpool is not open, computation occurs in serial
%                       mode. Default is 'default', meaning serial
%                       computation.
%    'UseSubstreams'  - Set to true to compute in parallel in a reproducible 
%                       fashion. Default is false. To compute reproducibly,
%                       set Streams to a type allowing substreams:
%                       'mlfg6331_64' or 'mrg32k3a'.
%          'Streams'  - These fields specify whether to perform clustering
%                       from multiple 'Start' values in parallel, and how
%                       to use random numbers when generating the starting
%                       points. For information on these fields see
%                       PARALLELSTATS.
%                       NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is FALSE,
%                       then the length of 'Streams' must equal the number of workers 
%                       used by KMEANS.  If a parallel pool is already open, this 
%                       will be the size of the parallel pool.  If a parallel pool 
%                       is not already open, then MATLAB may try to open a pool for 
%                       you (depending on your installation and preferences).
%                       To ensure more predictable results, it is best to use 
%                       the PARPOOL command and explicitly create a parallel pool 
%                       prior to invoking KMEANS with 'UseParallel' set to TRUE. 
%
%   'OnlinePhase' - Flag indicating whether KMEANS should perform an "on-line
%      update" phase in addition to a "batch update" phase.  The on-line phase
%      can be time consuming for large data sets, but guarantees a solution
%      that is a local minimum of the distance criterion, i.e., a partition of
%      the data where moving any single point to a different cluster increases
%      the total sum of distances.  'off' (the default) or 'on'.
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       opts = statset('Display','final');
%       [cidx, ctrs] = kmeans(X, 2, 'Distance','city', ...
%                             'Replicates',5, 'Options',opts);
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first phase
%   uses what the literature often describes as "batch" updates, where each
%   iteration consists of reassigning points to their nearest cluster
%   centroid, all at once, followed by recalculation of cluster centroids.
%   This phase occasionally (especially for small data sets) does not converge
%   to solution that is a local minimum, i.e., a partition of the data where
%   moving any single point to a different cluster increases the total sum of
%   distances.  Thus, the batch phase be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the second
%   phase.  The second phase uses what the literature often describes as
%   "on-line" updates, where points are individually reassigned if doing so
%   will reduce the sum of distances, and cluster centroids are recomputed
%   after each reassignment.  Each iteration during this second phase consists
%   of one pass though all the points.  The on-line phase will converge to a
%   local minimum, although there may be other local minima with lower total
%   sum of distances.  The problem of finding the global minimum can only be
%   solved in general by an exhaustive (or clever, or lucky) choice of
%   starting points, but using several replicates with random starting points
%   typically results in a solution that is a global minimum.
%
% References:
%
%   [1] Seber, G.A.F. (1984) Multivariate Observations, Wiley, New York.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York.

%   Copyright 1993-2014 The MathWorks, Inc.

if nargin < 2
    error(message('stats:kmeans:TooFewInputs'));
end

if ~isreal(X)
    error(message('stats:kmeans:ComplexData'));
end
wasnan = any(isnan(X),2);
hadNaNs = any(wasnan);
if hadNaNs
    warning(message('stats:kmeans:MissingDataRemoved'));
    X = X(~wasnan,:);
end

% n points in p dimensional space
[n, p] = size(X);

pnames = {   'distance'  'start' 'replicates' 'emptyaction' 'onlinephase' 'options' 'maxiter' 'display'};
dflts =  {'sqeuclidean' 'plus'          []  'singleton'         'off'        []        []        []};
[distance,start,reps,emptyact,online,options,maxit,display] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});

distNames = {'sqeuclidean','cityblock','cosine','correlation','hamming'};
distance = internal.stats.getParamVal(distance,distNames,'''Distance''');

switch distance
    case 'cosine'
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error(message('stats:kmeans:ZeroDataForCos'));
        end
        X =  bsxfun(@rdivide,X,Xnorm);
    case 'correlation'
        X = bsxfun(@minus, X, mean(X,2));
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error(message('stats:kmeans:ConstantDataForCorr'));
        end
        X =  bsxfun(@rdivide,X,Xnorm);
     case 'hamming'
       if  ~all( X(:) ==0 | X(:)==1)
            error(message('stats:kmeans:NonbinaryDataForHamm'));
      end
end

Xmins = [];
Xmaxs = [];
CC = [];
if ischar(start)
    startNames = {'uniform','sample','cluster','plus','kmeans++'};
    j = find(strncmpi(start,startNames,length(start)));
    if length(j) > 1
        error(message('stats:kmeans:AmbiguousStart', start));
    elseif isempty(j)
        error(message('stats:kmeans:UnknownStart', start));
    elseif isempty(k)
        error(message('stats:kmeans:MissingK'));
    end
    start = startNames{j};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error(message('stats:kmeans:UniformStartForHamm'));
        end
        Xmins = min(X,[],1);
        Xmaxs = max(X,[],1);
    end
elseif isnumeric(start)
    CC = start;
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1);
        error(message('stats:kmeans:StartBadRowSize'));
    elseif size(CC,2) ~= p
        error(message('stats:kmeans:StartBadColumnSize'));
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3);
        error(message('stats:kmeans:StartBadThirdDimSize'));
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
          CC = bsxfun(@minus, CC, mean(CC,2));
    end
else
    error(message('stats:kmeans:InvalidStart'));
end

emptyactNames = {'error','drop','singleton'};
emptyact = internal.stats.getParamVal(emptyact,emptyactNames,'''EmptyAction''');

[~,online] = internal.stats.getParamVal(online,{'on','off'},'''OnlinePhase''');
online = (online==1);

% 'maxiter' and 'display' are grandfathered as separate param name/value pairs
if ~isempty(display)
    options = statset(options,'Display',display);
end
if ~isempty(maxit)
    options = statset(options,'MaxIter',maxit);
end

options = statset(statset('kmeans'), options);
display = find(strncmpi(options.Display, {'off','notify','final','iter'},...
    length(options.Display))) - 1;
maxit = options.MaxIter;

if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k)==k))
    error(message('stats:kmeans:InvalidK'));
    % elseif k == 1
    % this special case works automatically
elseif n < k
    error(message('stats:kmeans:TooManyClusters'));
end

% Assume one replicate
if isempty(reps)
    reps = 1;
elseif ~internal.stats.isScalarInt(reps,1)
    error(message('stats:kmeans:BadReps'));
end

[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(options,true);

usePool = useParallel && poolsz>0;

% Prepare for in-progress
if display > 1 % 'iter' or 'final'
    if usePool
        % If we are running on a parallel pool, each worker will generate
        % a separate periodic report.  Before starting the loop, we
        % seed the parallel pool so that each worker will have an
        % identifying label (eg, index) for its report.
        internal.stats.parallel.distributeToPool( ...
            'workerID', num2cell(1:poolsz) );
        
        % Periodic reports behave differently in parallel than they do
        % in serial computation (which is the baseline).
        % We advise the user of the difference.
        
        if display == 3 % 'iter' only
            warning(message('stats:kmeans:displayParallel2'));
            fprintf('    worker\t  iter\t phase\t     num\t         sum\n' );
        end
    else
        if useParallel
            warning(message('stats:kmeans:displayParallel'));
        end
        if display == 3 % 'iter' only
            fprintf('  iter\t phase\t     num\t         sum\n');
        end
    end
end

if issparse(X) || ~isfloat(X) || strcmp(distance,'cityblock') || ...
        strcmp(distance,'hamming')
    [varargout{1:nargout}] = kmeans2(X,k, distance, emptyact,reps,start,...
        Xmins,Xmaxs,CC,online,display, maxit,useParallel, RNGscheme,usePool,...
        wasnan,hadNaNs,varargin{:});
    return;
end

emptyErrCnt = 0;

% Define the function that will perform one iteration of the
% loop inside smartFor
loopbody = @loopBody;

% Initialize nested variables so they will not appear to be functions here
totsumD = 0;
iter = 0;

X = X'; %Transpose data into column orientation
Xmins = Xmins';
Xmaxs = Xmaxs';

% Perform KMEANS replicates on separate workers.
ClusterBest = internal.stats.parallel.smartForReduce(...
    reps, loopbody, useParallel, RNGscheme, 'argmin');

% Extract the best solution
varargout{1} = ClusterBest{5};%idxbest = ClusterBest{5};
varargout{2} = ClusterBest{6}';%Cbest = ClusterBest{6};
varargout{3} = ClusterBest{3}; %sumDbest = ClusterBest{3};
totsumDbest = ClusterBest{1};

if nargout > 3
    varargout{4} = ClusterBest{7}; %Dbest
end

if display > 1 % 'final' or 'iter'
    fprintf('%s\n',getString(message('stats:kmeans:FinalSumOfDistances',sprintf('%g',totsumDbest))));
end
 
if hadNaNs
    varargout{1} = statinsertnan(wasnan, varargout{1});% idxbest 
    if nargout > 3
        varargout{4} = statinsertnan(wasnan, varargout{4}); %Dbest
    end
end

    function cellout = loopBody(rep,S)
        
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        
        if display > 1 % 'iter'
            if usePool
                dispfmt = '%8d\t%6d\t%6d\t%8d\t%12g\n';
                labindx = internal.stats.parallel.workerGetValue('workerID');
            else
                dispfmt = '%6d\t%6d\t%8d\t%12g\n';
            end
        end
        
        cellout = cell(7,1);  % cellout{1} = total sum of distances
                              % cellout{2} = replicate number
                              % cellout{3} = sum of distance for each cluster
                              % cellout{4} = iteration
                              % cellout{5} = idx;
                              % cellout{6} = Center
                              % cellout{7} = Distance
        
        % Populating total sum of distances to Inf. This is used in the
        % reduce operation if update fails due to empty cluster.
        cellout{1} = Inf;
        cellout{2} = rep;
        
        switch start
            case 'uniform'
                %C = Xmins(:,ones(1,k)) + rand(S,[p,k]).*(Xmaxs(:,ones(1,k))-Xmins(:,ones(1,k)));
                C = Xmins(:,ones(1,k)) + rand(S,[k,p])'.*(Xmaxs(:,ones(1,k))-Xmins(:,ones(1,k)));
                % For 'cosine' and 'correlation', these are uniform inside a subset
                % of the unit hypersphere.  Still need to center them for
                % 'correlation'.  (Re)normalization for 'cosine'/'correlation' is
                % done at each iteration.
                if isequal(distance, 'correlation')
                    C = bsxfun(@minus, C, mean(C,1));
                end
                if isa(X,'single')
                    C = single(C);
                end
            case 'sample'
                C = X(:,randsample(S,n,k));
            case 'cluster'
                Xsubset = X(:,randsample(S,n,floor(.1*n)));
                % Turn display off for the initialization
                optIndex = find(strcmpi('options',varargin));
                if isempty(optIndex)
                    opts = statset('Display','off');
                    varargin = [varargin,'options',opts];
                else
                    varargin{optIndex+1}.Display = 'off';
                end
                [~, C] = kmeans(Xsubset', k, varargin{:}, 'start','sample', 'replicates',1);
                C = C';
            case 'numeric'
                C = CC(:,:,rep)';
                if isa(X,'single')
                    C = single(C);
                end
            case {'plus','kmeans++'}
                % Select the first seed by sampling uniformly at random
                index = zeros(1,k);
                [C(:,1), index(1)] = datasample(S,X,1,2);
                minDist = inf(n,1);
           
                % Select the rest of the seeds by a probabilistic model
               for ii = 2:k                    
                    minDist = min(minDist,distfun(X,C(:,ii-1),distance));
                    denominator = sum(minDist);
                    if denominator==0 || isinf(denominator) || isnan(denominator)
                        C(:,ii:k) = datasample(S,X,k-ii+1,2,'Replace',false);
                        break;
                    end
                    sampleProbability = minDist/denominator;
                    [C(:,ii), index(ii)] = datasample(S,X,1,2,'Replace',false,...
                        'Weights',sampleProbability);        
                end
        end
        if ~isfloat(C)      % X may be logical
            C = double(C);
        end
        
        % Compute the distance from every point to each cluster centroid and the
        % initial assignment of points to clusters
        D = distfun(X, C, distance, 0, rep, reps);
        [d, idx] = uniform_min(D);
        %[d,idx] = min(D,[],2);
        m = accumarray(idx,1,[k,1])';
        
        try % catch empty cluster errors and move on to next rep
            
            % Begin phase one:  batch reassignments
            converged = batchUpdate();
            
            % Begin phase two:  single reassignments
            if online
                converged = onlineUpdate();
            end
            
            
            if display == 2 % 'final'
                fprintf('%s\n',getString(message('stats:kmeans:IterationsSumOfDistances',rep,iter,sprintf('%g',totsumD) )));
            end
            
            if ~converged
                if reps==1
                    warning(message('stats:kmeans:FailedToConverge', maxit));
                else
                    warning(message('stats:kmeans:FailedToConvergeRep', maxit, rep));
                end
            end
            
            % Calculate cluster-wise sums of distances
            nonempties = find(m>0);
            D(:,nonempties) = distfun(X, C(:,nonempties), distance, iter, rep, reps);
            d = D((idx-1)*n + (1:n)');
            sumD = accumarray(idx,d,[k,1]);
            totsumD = sum(sumD(nonempties));
            
            % Save the best solution so far
            cellout = {totsumD,rep,sumD,iter,idx,C,D}';
           
            % If an empty cluster error occurred in one of multiple replicates, catch
            % it, warn, and move on to next replicate.  Error only when all replicates
            % fail.  Rethrow an other kind of error.
        catch ME
            if reps == 1 || (~isequal(ME.identifier,'stats:kmeans:EmptyCluster')  && ...
                         ~isequal(ME.identifier,'stats:kmeans:EmptyClusterRep'))
                rethrow(ME);
            else
                emptyErrCnt = emptyErrCnt + 1;
                warning(message('stats:kmeans:EmptyClusterInBatchUpdate', rep, iter));
                if emptyErrCnt == reps
                    error(message('stats:kmeans:EmptyClusterAllReps'));
                end
            end
        end % catch
        
        %------------------------------------------------------------------
        
        function converged = batchUpdate()
            
            % Every point moved, every cluster will need an update
            moved = 1:n;
            changed = 1:k;
            previdx = zeros(n,1);
            prevtotsumD = Inf;
            
            %
            % Begin phase one:  batch reassignments
            %
            
            iter = 0;
            converged = false;
            while true
                iter = iter + 1;
                
                % Calculate the new cluster centroids and counts, and update the
                % distance from every point to those new cluster centroids
                [C(:,changed), m(changed)] = gcentroids(X, idx, changed, distance);
                D(:,changed) = distfun(X, C(:,changed), distance, iter, rep, reps);
                
                % Deal with clusters that have just lost all their members
                empties = changed(m(changed) == 0);
                if ~isempty(empties)
                    if strcmp(emptyact,'error')
                        if reps==1
                            error(message('stats:kmeans:EmptyCluster', iter));
                        else
                            error(message('stats:kmeans:EmptyClusterRep', iter, rep));
                        end
                    end
                    switch emptyact
                        case 'drop'
                            if reps==1
                                warning(message('stats:kmeans:EmptyCluster', iter));
                            else
                                warning(message('stats:kmeans:EmptyClusterRep', iter, rep));
                            end
                            % Remove the empty cluster from any further processing
                            D(:,empties) = NaN;
                            changed = changed(m(changed) > 0);
                        case 'singleton'
                            for i = empties
                                d = D((idx-1)*n + (1:n)'); % use newly updated distances
                                
                                % Find the point furthest away from its current cluster.
                                % Take that point out of its cluster and use it to create
                                % a new singleton cluster to replace the empty one.
                                [~, lonely] = max(d);
                                from = idx(lonely); % taking from this cluster
                                if m(from) < 2
                                    % In the very unusual event that the cluster had only
                                    % one member, pick any other non-singleton point.
                                    from = find(m>1,1,'first');
                                    lonely = find(idx==from,1,'first');
                                end
                                C(:,i) = X(:,lonely);
                                m(i) = 1;
                                idx(lonely) = i;
                                D(:,i) = distfun(X, C(:,i), distance, iter, rep, reps);
                                
                                % Update clusters from which points are taken
                                [C(:,from), m(from)] = gcentroids(X, idx, from, distance);
                                D(:,from) = distfun(X, C(:,from), distance, iter, rep, reps);
                                changed = unique([changed from]);
                            end
                    end
                end
                
                % Compute the total sum of distances for the current configuration.
                totsumD = sum(D((idx-1)*n + (1:n)'));
                % Test for a cycle: if objective is not decreased, back out
                % the last step and move on to the single update phase
                if prevtotsumD <= totsumD
                    idx = previdx;
                    [C(:,changed), m(changed)] = gcentroids(X, idx, changed, distance);
                    iter = iter - 1;
                    break;
                end
                if display > 2 % 'iter'
                    if usePool
                        fprintf(dispfmt,labindx,iter,1,length(moved),totsumD);
                    else
                        fprintf(dispfmt,iter,1,length(moved),totsumD);
                    end
                end
                if iter >= maxit
                    break;
                end
                
                % Determine closest cluster for each point and reassign points to clusters
                previdx = idx;
                prevtotsumD = totsumD;
                [d, nidx] = uniform_min(D);
                %[d, nidx] = min(D,[],2);
                
                % Determine which points moved
                moved = find(nidx ~= previdx);
                if ~isempty(moved)
                    % Resolve ties in favor of not moving
                    moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
                end
                if isempty(moved)
                    converged = true;
                    break;
                end
                idx(moved) = nidx(moved);
                
                % Find clusters that gained or lost members
                changed = unique([idx(moved); previdx(moved)])';
                
            end % phase one
            
        end % nested function
        
        %------------------------------------------------------------------
        
        function converged = onlineUpdate()
                       
            %
            % Begin phase two:  single reassignments
            %
            changed = find(m > 0);
            lastmoved = 0;
            nummoved = 0;
            iter1 = iter;
            converged = false;
            Del = NaN(n,k); % reassignment criterion
            while iter < maxit
                % Calculate distances to each cluster from each point, and the
                % potential change in total sum of errors for adding or removing
                % each point from each cluster.  Clusters that have not changed
                % membership need not be updated.
                %
                % Singleton clusters are a special case for the sum of dists
                % calculation.  Removing their only point is never best, so the
                % reassignment criterion had better guarantee that a singleton
                % point will stay in its own cluster.  Happily, we get
                % Del(i,idx(i)) == 0 automatically for them.
                switch distance
                    case 'sqeuclidean'
                        for i = changed
                            mbrs = (idx == i)';
                            sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                            if m(i) == 1
                                sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
                            end
                          Del(:,i) = (m(i) ./ (m(i) + sgn)) .* sum((bsxfun(@minus, X, C(:,i))).^2, 1);
                        end
                      case {'cosine','correlation'}
                        % The points are normalized, centroids are not, so normalize them
                        normC = sqrt(sum(C.^2, 1));
                        if any(normC < eps(class(normC))) % small relative to unit-length data points
                            if reps==1
                                error(message('stats:kmeans:ZeroCentroid', iter));
                            else
                                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
                            end
                            
                        end
                        % This can be done without a loop, but the loop saves memory allocations
                        for i = changed
                            XCi =  C(:,i)'*X;
                            mbrs = (idx == i)';
                            sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
                            Del(:,i) = 1 + sgn .*...
                                (m(i).*normC(i) - sqrt((m(i).*normC(i)).^2 + 2.*sgn.*m(i).*XCi + 1));
                        end
                end
                
                % Determine best possible move, if any, for each point.  Next we
                % will pick one from those that actually did move.
                previdx = idx;
                prevtotsumD = totsumD;
                [minDel, nidx] = min(Del, [], 2);
                moved = find(previdx ~= nidx);
                moved(m(previdx(moved))==1)=[];
                if ~isempty(moved)
                    % Resolve ties in favor of not moving
                    moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
                end
                if isempty(moved)
                    % Count an iteration if phase 2 did nothing at all, or if we're
                    % in the middle of a pass through all the points
                    if (iter == iter1) || nummoved > 0
                        iter = iter + 1;
                        if display > 2 % 'iter'
                            if usePool
                                fprintf(dispfmt,labindx,iter,2,length(moved),totsumD);
                            else
                                fprintf(dispfmt,iter,2,length(moved),totsumD);
                            end
                        end
                    end
                    converged = true;
                    break;
                end
                
                % Pick the next move in cyclic order
                moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
                
                % If we've gone once through all the points, that's an iteration
                if moved <= lastmoved
                    iter = iter + 1;
                    if display > 2 % 'iter'
                        if usePool
                            fprintf(dispfmt,labindx,iter,2,length(moved),totsumD);
                        else
                            fprintf(dispfmt,iter,2,length(moved),totsumD);
                        end
                    end
                    if iter >= maxit, break; end
                    nummoved = 0;
                end
                nummoved = nummoved + 1;
                lastmoved = moved;
                
                oidx = idx(moved);
                nidx = nidx(moved);
                totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
                
                % Update the cluster index vector, and the old and new cluster
                % counts and centroids
                idx(moved) = nidx;
                m(nidx) = m(nidx) + 1;
                m(oidx) = m(oidx) - 1;
                switch distance
                    case {'sqeuclidean','cosine','correlation'}
                        C(:,nidx) = C(:,nidx) + (X(:,moved) - C(:,nidx)) / m(nidx);
                        C(:,oidx) = C(:,oidx) - (X(:,moved) - C(:,oidx)) / m(oidx);
                end
                changed = sort([oidx nidx]);
            end % phase two
            
        end % nested function
        
    end

end % main function

%------------------------------------------------------------------

function D = distfun(X, C, dist, iter,rep, reps)
%DISTFUN Calculate point to cluster centroid distances.

switch dist
    case 'sqeuclidean'
        D = pdist2(X',C','squaredeuclidean');  
    case {'cosine','correlation'}
        % The points are normalized, centroids are not, so normalize them
        normC = sqrt(sum(C.^2, 1));
        if any(normC < eps(class(normC))) % small relative to unit-length data points
            if reps==1
                error(message('stats:kmeans:ZeroCentroid', iter));
            else
                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
            end
            
        end
        C = bsxfun(@rdivide,C,normC);
        D = pdist2(X',C','cosine');  
end
end % function

%------------------------------------------------------------------
function [centroids, counts] = gcentroids(X, index, clusts, dist)
%GCENTROIDS Centroids and counts stratified by group.
p = size(X,1);
num = length(clusts);

centroids = NaN(p,num,'like',X);
counts = zeros(1,num,'like',X);

for i = 1:num
    members = (index == clusts(i));
    if any(members)
       counts(i) = sum(members);
       switch dist
           case {'sqeuclidean','cosine','correlation'}
              centroids(:,i) = sum(X(:,members),2) / counts(i);
      end
    end
end
end % function

function [d, idx] = uniform_min(D); % cout modification
max_size = floor(size(D,1)/size(D,2));
idx = zeros(size(D,1),1);
d = idx;
clustn = zeros(1,size(D,2));
j=1;
while j<size(D,1)+1
    idx2 = idx;
    use_cluster = find(clustn<=max_size);
    id = idx(~idx);
    dd = d(~idx);
    
    [d2,i2] = min(D(~idx,use_cluster),[],2);
    i2 = use_cluster(i2);
    [d1,i1] = min(d2,[],1);
    id(i1) = i2(i1);
    dd(i1(1)) = d1;
    
    idx(~idx2) = id;
    d(~idx2) = dd;
    
    clustn(i2(i1)) = clustn(i2(i1))+1;
    j = j+1;
end
end
