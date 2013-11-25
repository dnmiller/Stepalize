function [A, B, C, D] = stepalize(y, varargin)
% stepalize.m: Build a state-space realization from a step-response.
% 
% This uses a step-based realization (SBR) algorithm to create a linear,
% discrete-time, state-space system from a (possibly multivariable)
% step-response measurement.
% 
%   [A, B, C, D] = stepalize(y, args) 
% 
%       Create a step-response model from a step-response measurement y. If
%       y is a matrix, than each column of y is treated as a separate
%       output signal. 
% 
%       args are optional 'Property', 'Value' pairs of the following:
%       
%       'Order' - Order of the model 
%           - If not specified, the user will be presented with a
%             singular-value plot on which they may select the system
%             order.
% 
%       'BlockRows' - Number of block rows for data matrices
%           - If not specified, then the number is either
%               (1) min(floor(size(y, 1)/4), 50) if the order is not
%                   specified
%               (2) n * 4 where n is the specified order
% 
%       'TimeDelay' - Number of time delays in data
%           - By default, 1 time delay is assumed, and D is always 0. 
%           - Set to 0 to estimate a nonzero D.
%           - If > 1, then the data is shifted so that there is 1 time
%             delay. D is again 0.
%
%       'Eigenvalues' - Constraints to place on eigenvalues
%           This has the following possible values:
% 
%           - 'Stable' - force stable eigenvalues
%               - Requires all eigenvalues to satisfy |z| < 1.0 - ds.
%               - ds may be optionally set by the 'DeltaS' property.
%               - By default, ds = 1e-4.
% 
%           - 'Real' - force strictly real eigenvalues
%               - Requires all eigenvalues to satisfy |Im(z)| < dr.
%               - dr may be optionally set by the 'DeltaR' property.
%               - By default, dr = 1e-8.
% 
%           - 'Positive' - force eigenvalues with positive real parts
%               - Requires all eigenvalues to satisfy Re(z) > dp.
%               - dp may be optionally set by the 'DeltaP' property.
%               - By default, dp = 1e-4.
% 
%           If one is used, it may be passed as a single string, i.e.
%               
%               ..., 'Eigenvalues', 'Real', ...
% 
%           For multiple choices, a cell array may be used, i.e.
% 
%               ..., 'Eigenvalues', {'Real', 'Stable'}, ...
% 
%       'SteadyState' - Guarantee the model has a fixed stead-state value.
%           The value for this property should be a 1 x n_y vector which
%           will be settling value for the model's step response.
% 
%       'NoOvershoot' - Guarantee the model does not overshoot the
%           steady-state value of the data.
%           - This is often infeasible, depending on the eigenvalues. If a
%             step-response with no overshoot is absolutely required,
%             constrain the eigenvalues to be real.
% 
%       'NoUndershoot' - Restrict resulting model to have no undershoot in
%           its step response. 
%           - This generally eliminates non-minimum-phase behavior from the
%             model.
%           - Undershoot is defined as when the step response has the
%             opposite sign of the steady-state value for a given output
%             channel. This allows for the step response to still have
%             values < 0, but only when the steady-state value at that
%             channel is < 0 as well.
% 
%        - NOTE: Combinations of 'SteadyState', 'NoOvershoot', and
%          'NoUndershoot' are prone to overconstraining the model and will
%          often result in an infeasibility error or the solution B = 0.
%          Only use them in combination when the estimate is slightly off.
%          Feasibility issues are also reduced when the system is
%          parameterized with a feed-through (D-matrix) term.
% 
%       'UseAC' - Use pre-calculated matrices A and C to compute B and D
%           only.
%           - The value for this must be a 2-element cell array with the
%             first element equal to A and the second equal to C. 
%           - Options for calculating B and D are still available.
%           - Any options affecting the computation of A and C are ignored.
% 
%       'WeightBD' - Weight the solution for B and D. 
%           - This should be a N x 1 vector for which element i corresponds
%             to the weight of the solution for sample i.
%           - This is not available if constraints are used with B and D. A
%             warning will be generated.
%       
% Example:
% 
%   % Generate an n-th order model with a feed-through term.
%   [A, B, C, D] = stepalize(y, 'Order', n, 'TimeDelay', 0)
% 
% References:
% 
%   D. Miller and R.A. de Callafon, "Step-Based Realization with Eigenvalue
%   Constraints," 16th IFAC Symposium on System Identification, Brussels:
%   2012.


%   (C) 2012 D. Miller, UCSD.
    args = parse_args(y, varargin{:});
    
    % If the user requests verbose output, then print things to the
    % command window as we progress.
    if args.Verbose
        log_fnc = @(x) fprintf(x);
    else
        log_fnc = @(x) x;
    end
    
    % We just shift for time delays.
    if args.TimeDelay > 1
        args.y = args.y(args.TimeDelay:end, :);
    end

    % Estimate the system.
    if isempty(args.UseAC)
        [A, C] = solveAC(args, log_fnc);
    else
        A = args.UseAC{1};
        C = args.UseAC{2};
    end
    [B, D] = solveBD(args, A, C, log_fnc);
end


function [A, C] = solveAC(args, log_fnc)
    % Output signal dimension.
    ny = size(args.y, 2);
    
    % Construct the output data matrices.
    log_fnc('Constructing data matrices....'); tic;
    
    % Make up a default block row size if none was given. No good way to do
    % this.
    if isempty(args.BlockRows)
        if isempty(args.Order)
            r = min(floor(size(args.y, 1)/4), 50);
        else
            r = floor(args.Order*4);
        end
    else
        r = args.BlockRows;
    end
    
    Y = datahankel(args.y(2:end, :), r);
    Y0 = Y(1:end-ny, :);
    Y1 = Y(ny+1:end, :);
    
    % Construct the column-identical input-product matrices.
    log_fnc('....');
    yx = args.y(1:r, :)';
    yx = yx(:);
    M = repmat(yx(:), 1, size(Y, 2));
    M0 = M(1:end-ny, :);
    M1 = M(ny+1:end, :);

    % Construct the weighted Hankel matrices.
    log_fnc('....');
    R0 = Y0 - M0;
    R1 = Y1 - M1;

    clear Y Y0 Y1 M M0 M1;
    log_fnc(['...Elapsed time ' num2str(toc), 's\n']);
    
    % Take the SVD of the approximate weighted Hankel matrix.
    log_fnc('Taking SVD...'); tic;
    [Ux, Sx, Vx] = svd(R0);
    log_fnc(['............................Elapsed time ' num2str(toc), 's\n']);
    
    % Choose the order if necessary.
    if isempty(args.Order)
        n = selectorder(diag(Sx));
        if n < 1
            ex = MException('STEPALIZE:BadOrderSelect', ...
                    'Invalid order selection. Try again.');
            throw(ex);
        end
    else
        n = args.Order;
    end
    Un = Ux(:, 1:n);
    Sn = Sx(1:n, 1:n);
    Vn = Vx(:, 1:n);

    % Estimate observability / controllability matrices.
    Ob = Un*Sn^(1/2);
    Ob_i = Sn^(-1/2)*Un';
    CrUp_i = Vn*Sn^(-1/2);

    % C is taken from the first ny rows of the estimate of the
    % extended observability matrix.
    C = Ob(1:ny, :);

    % Compute the state dynamics (A matrix).
    log_fnc('Computing state dynamics...'); tic;
    
    A = Ob_i*R1*CrUp_i;
    
    if ~isempty(args.Eigenvalues) && ~satisfies_eig_constraints(A, args)
        % If there are constraints on the eigenvalues, we must solve an
        % SDP.
        [P, Q, M] = build_LMIs(n, args);
        assign(Q, A);
        
        Constraints = [P > 0; M > 0];
        Costs = norm(Q - Ob_i*R1*CrUp_i*P, 'fro');
        
        log_fnc('\n');
        
        options = sdpsettings(...
            'solver', 'sdpt3', ...
            'verbose', args.Verbose, ...
            'usex0', 1);
        
        d = solvesdp(Constraints, Costs, options);
        
        if d.problem ~= 0
            ex = MException('STEPALIZE:YALMIPError', ...
                ['YALMIP reported error: "' d.info '" when solving for A.']);
            throw(ex);
        end
        
        % Some cases, when testing, we want to make sure that the estimate
        % is consistent for the deterministic case. The SDP however will
        % sometimes provide an A with the correct eigenvalues but in a
        % basis different than that of C. This changes C to the correct
        % basis.
        %
        % Normally, however, we don't enter the optimization routine when
        % the eigenvalues of the unconstrained estimate already satisfy the
        % constraints. If the eigenvalues of the two estimates are
        % different, then we don't want to change the basis of C, since
        % this will introduce some errors.
        %
        % This should only be uncommented for testing purposes.
        %
%         [S, ~] = eig(double(P)\double(Q));
%         [Su, ~] = eig(A);
%         
%         T = real(S/Su);
%         C = C/T;

        A = double(P)\double(Q);
    end
    log_fnc(['..............Elapsed time ' num2str(toc), 's\n']);
end


function [B, D] = solveBD(args, A, C, log_fnc)
    [ny, n] = size(C);
    nu = 1;
    N = size(args.y, 1);
    
    % Can't caluclate B and D for unstable A. We don't raise an exception
    % here because A and C can still be theoretically correct.
    if any(abs(eig(A)) >= 1)
        warning('STEPALIZE:unstableA', ...
            'Model unstable. B and D could not be identified');
        B = NaN*ones(n, nu);
        D = NaN*ones(ny, nu);
        return;
    end
    
    % Compute the regressors.
    [phiB, phiD] = calcPhi(args, A, C, n, N, ny, nu, log_fnc);
        
    % Find the least-squares solution for B, D.
    if args.NoOvershoot || args.NoUndershoot || ~isempty(args.SteadyState)
        phiSS = C/(eye(n) - A);
        [B, D] = solveBD_con(args, phiB, phiD, phiSS, N, n, ny, nu);
    else
        phi = [phiB', phiD'];
        [B, D] = solveBD_uncon(args, phi, n, ny, nu, log_fnc);
    end
end


function [B, D] = solveBD_uncon(args, phi, n, ny, nu, log_fnc)
% Solve the regression for B and D in the unconstrained case.

    % Vectorize y to use in the least-squares problem.
    y = args.y';
    y = y(:);

    log_fnc('Computing the unconstrained regressor...'); tic;
    if isempty(args.WeightBD)
        theta = phi\y;
    else
        theta = lscov(phi, y, args.WeightBD);
    end

    B = reshape(theta(1:n*nu), n, nu);

    if args.TimeDelay == 0
        D = reshape(theta(n*nu+1:end), ny, nu);
    else
        D = zeros(ny, nu);
    end
    log_fnc(['..............Elapsed time ' num2str(toc), 's\n']);
end


function [B, D] = solveBD_con(args, phiB, phiD, phiSS, N, n, ny, nu)
% Solve the regression for B and D in the constrained case.
% 
% The constraints naturally make this a little more complicated. Since the
% total regressor must be built up using sdpvar objects, each regressor for
% B and D must be passed in separately.
% 
% Additionally, if the steady-state value is constrained, then C(I - A)^-1
% must also be passed in. This is phiSS in the function arguments.

    % This is how we build up the total regressor. We need to keep phi
    % linked to B and D for constraint purposes later on.
    B = sdpvar(n, nu, 'full');
    ye = phiB'*B;

    if args.TimeDelay == 0
        D = sdpvar(ny, nu, 'full');
        ye = ye + phiD'*D;
    else
        D = zeros(ny, nu);
    end

    % Vectorize y to use in the least-squares problem.
    y = args.y';
    y = y(:);

    % Determine the sign of the steady-state value of each output. We'll
    % use this to enforce response constraints.
    ny_sign = zeros(ny, 1);
    if isempty(args.SteadyState)
        for k = 1:ny
            ny_sign(k) = sign(args.y(end, k));
        end
    else
        for k = 1:ny
            ny_sign(k) = sign(args.SteadyState(k));
        end            
    end
    
    % Normalize the signs of output signals so that everything is positive.
    % This results in "undershoot" and "overshoot" being measured relative
    % to the steady-state instead of relative to 0.
    sign_norm = repmat(ny_sign, N, 1);
    y = y.*sign_norm;
    ye = ye.*sign_norm;

    % There are various constraints that can be put on the estimate, such
    % as no undershoot, no overshoot, and fixed steady-state value.
    Constraints = [];
    
    % No undershoot is simple enough. We've already changed the signs in
    % ye, so we can just make it > 0.
    if args.NoUndershoot
        Constraints = [Constraints, ye > 0];
    end

    % No overshoot involves the steady-state, which is harder.
    if args.NoOvershoot
        % No good way to compute the steady-state here automatically that
        % doesn't overconstrain the system. So we either look for user
        % input or use the final measurement of the data for 
        if isempty(args.SteadyState)
            yss = args.y(end, :)';
        else
            yss = args.SteadyState';
        end
        Constraints = [Constraints, ye < repmat(yss, size(ye, 1)/ny, 1)];
    end

    % Add steady-state constraints. 
    if ~isempty(args.SteadyState)
        Constraints = [Constraints, args.SteadyState' == phiSS*B + D];
    end

    % Trying to weight the constrained case causes numerical problems. This
    % may be a limitation of the convex optimization solvers, and there is
    % probably nothing we can do about it.
    if ~isempty(args.WeightBD)
        warning('on', 'STEPALIZE:NoWeighting');
        warning('STEPALIZE:NoWeighting', ...
            'Weighting of BD solution currently unavailable with constraints.');
    end

    Cost = norm(y - ye);

    options = sdpsettings(...
        'solver', 'sdpt3', ...
        'verbose', args.Verbose);

    d = solvesdp(Constraints, Cost, options);
    if d.problem ~= 0
        ex = MException('STEPALIZE:YALMIPError', ...
            ['YALMIP reported error: "' d.info '" when solving for B and D.']);
        throw(ex);
    end

    B = double(B);
    D = double(D);
end

function [phiB, phiD] = calcPhi(args, A, C, n, N, ny, nu, log_fnc)
% This calculates the regressors to solve the linear least-squares problem
% that finds B and D. This regressor is used in both the constrained
% and unconstrained versions.

    % The regressor for B is actually a state sequence of a dual system,
    % but with each input channel considered independently. See ref [1] for
    % details.
    log_fnc('Constructing input dynamics regressor....'); tic;
    phiB = zeros(n*nu, ny*N);
    for i = 1:nu
        for j = 1:ny
            uk = [zeros(N, j-1), ones(N, 1), zeros(N, ny-j)];
            x = ltitr(A', C', uk);
            
            % Copy the state sequence into the regressor for B.
            for k = 1:N
                phiB(n*(i-1)+1:n*i, ny*(k-1)+j) = x(k, :)';
            end
        end
    end
    log_fnc(['Elapsed time ' num2str(toc), 's\n']);
    
    % Calculate the regressor for D if needed.
    u = ones(N, 1);
    if args.TimeDelay == 0
        log_fnc('Constructing feed-through regressor...'); tic;
        
        phiD = zeros(nu*ny, ny*N);
        for k = 1:N
            % This for-loop takes advantage of the fact that we're taking a
            % kroneker product with an identity, so the matrix has a
            % special structure that lets us avoid use of the Matlab kron
            % function, which is horribly slow.
            for i = 1:nu
                phiD(ny*(i-1)+1:ny*i, ny*(k-1)+1:ny*k) = eye(ny)*u(k, i);
            end
        end
        log_fnc(['..............Elapsed time ' num2str(toc), 's\n']);
    else
        phiD = [];
    end
end


function r = satisfies_eig_constraints(A, args)
% Given a matrix A, determine if it is in the LMI region described by the
% constraints.

    if isempty(args.DeltaS)
        ds = 1e-4;
    else
        ds = args.DeltaS;
    end
    
    if isempty(args.DeltaR)
        dr = 1e-8;
    else
        dr = args.DeltaR;
    end
    
    if isempty(args.DeltaP)
        dp = 1e-4;
    else
        dp = args.DeltaP;
    end
    
    % These describe "LMI regions," which mean that if a complex number is
    % in some convex region of the complex plane, then f(z) > 0 (f(z) is a
    % matrix function).
    
    if any(strcmpi(args.Eigenvalues, 'Stable'));
        alpha_s = (1 - ds)*eye(2);
        beta_s = [0, 1; 0, 0];
        f_s = @(z) alpha_s + beta_s*z + beta_s'*z';
    else
        f_s = @(z) [];
    end

    if any(strcmpi(args.Eigenvalues, 'Real'))
        alpha_r = dr*eye(2);
        beta_r = [0, 0.5; -0.5, 0];
        f_r = @(z) alpha_r + beta_r*z + beta_r'*z';
    else
        f_r = @(z) [];
    end

    if any(strcmpi(args.Eigenvalues, 'Positive'));
        alpha_p = dp*[2, 0; 0, -2];
        beta_p = [0, 0; 0, 1];
        f_p = @(z) alpha_p + beta_p*z + beta_p'*z';
    else
        f_p = @(z) [];
    end
    
    f = @(z) blkdiag(f_s(z), f_r(z), f_p(z));

    lambda = eig(A);
    r = true;
    for k = 1:length(lambda)
        if any(eig(f(lambda(k))) < 0)
            r = false;
            break;
        end
    end
end


function [P, Q, M] = build_LMIs(n, args)
% Build the constraint matrix that describes the convex region the
% eigenvalues are constrained to.

    if isempty(args.DeltaS)
        ds = 1e-4;
    else
        ds = args.DeltaS;
    end
    
    if isempty(args.DeltaR)
        dr = 1e-8;
    else
        dr = args.DeltaR;
    end
    
    if isempty(args.DeltaP)
        dp = 1e-4;
    else
        dp = args.DeltaP;
    end
    
    P = sdpvar(n);
    Q = sdpvar(n, n, 'full');

    % This is our constraint matrix. It is a block-diagonal matrix
    % in which each constraint corresponds to one block.
    M = [];

    if any(strcmpi(args.Eigenvalues, 'Stable'));
        M = blkdiag(M, [
                (1 - ds)*P,     Q;
                Q',             (1 - ds)*P]);
    end

    if any(strcmpi(args.Eigenvalues, 'Real'))
        M = blkdiag(M, [
                dr*P,           0.5*(Q' - Q)
                0.5*(Q - Q'),   dr*P]);
    end

    if any(strcmpi(args.Eigenvalues, 'Positive'));
        M = blkdiag(M, [
                2*dp*P,         zeros(n);
                zeros(n),       Q + Q' - 2*dp*P]);
    end
end


function args = parse_args(y, varargin)
% Parse function arguments.
    
    % We'll use Matlab's built-in inputParser object to parse the input
    % arguments. These are returned in a giant structure to the calling
    % function.
    isnumericmat = @(x) validateattributes(x, ...
        {'numeric'}, {'real', 'nonempty', 'nonnan', 'finite'});
    
    isdimension = @(x) validateattributes(x, ...
        {'numeric'}, {'positive', 'real', 'scalar', 'integer'});

    isrealvector = @(x) validateattributes(x, ...
        {'numeric'}, {'real', 'nonnan', 'finite', 'vector'});
    
    isdelta = @(x) validateattributes(x, ...
        {'numeric'}, {'real', 'scalar', 'finite', 'nonnan', 'positive'});
    
    istimedelay = @(x) validateattributes(x, ...
        {'numeric'}, {'real', 'nonnegative', 'integer', 'scalar'});

    p = inputParser;
    p.addRequired('y', isnumericmat)
    p.addParamValue('Order', [], isdimension);
    p.addParamValue('BlockRows', [], isdimension);
    p.addParamValue('Eigenvalues', [], @validate_eig_constraints);
    p.addParamValue('DeltaP', [], isdelta);
    p.addParamValue('DeltaR', [], isdelta);
    p.addParamValue('DeltaS', [], isdelta);
    p.addParamValue('TimeDelay', 1, istimedelay);
    p.addParamValue('SteadyState', [], isrealvector);
    p.addParamValue('NoUndershoot', false, @(x) islogical(x));
    p.addParamValue('NoOvershoot', false, @(x) islogical(x));
    p.addParamValue('Verbose', false, @(x) islogical(x));
    p.addParamValue('WeightBD', [], isrealvector);
    p.addParamValue('UseAC', {}, @(x) iscell(x));

    p.parse(y, varargin{:});

    args = p.Results;
    
    % We know now that y is valid, so check steady-state value.
    if ~isempty(args.SteadyState)
        try
            validateattributes(args.SteadyState, {'numeric'}, ...
                                             {'size', [1, size(args.y, 2)]});
        catch ex
            if strcmp(ex.identifier, 'MATLAB:incorrectSize')
                err = MException('STEPALIZE:incorrectSteadyStateSize', ...
                    'SteadyState is incorrect size.');
                throw(err);
            end
        end
    end
    
%     if ~isempty(args.WeightBD)
%         args.WeightBD
    
    % Force the Eigenvalues field to be a cell array, even if the user
    % entered a string.
    if ischar(args.Eigenvalues)
        args.Eigenvalues = {args.Eigenvalues};
    end
    
    % YALMIP is required for eigenvalue constraints.
    if ~isempty(args.Eigenvalues) && exist('yalmip', 'dir') ~= 7
        err = MException('STEPALIZE:noYALMIP', ...
            'YALMIP is required for eigenvalue constraints.');
        err.throw();
    end
end


function res = validate_eig_constraints(c)
% Validate the arguments for eigenvalue constraints. This is meant to be
% passed as a function handle to Matlab's inputParser object.

    valid = {'', 'Real', 'Positive', 'Stable'};
    reqstr = 'Constraints must be a string or a 1-D cell of strings.';
    
    res = true;
    % If the constraint is a character array, check to see if it has a
    % valid value.
    if ischar(c)
        if any(strcmpi(c, valid))
            res = true;
        else
            disp([c ' is not a valid eigenvalue constraint.']);
            res = false;
        end
        
    % If the constraint is a cell array, iterate over each, making sure
    % the entries have valid values.
    elseif iscell(c) && isvector(c)
        for k = 1:length(c)
            if ~ischar(c{k})
                disp(reqstr);
                res = false;
            elseif ~any(strcmpi(c{k}, valid))
                disp([c{k} ' is not a valid eigenvalue constraint.']);
                res = false;
            end
        end
    
    % Otherwise fail.
    else
        disp(reqstr);
        res = false;
    end
end


% -------------------------------------------------------------------------
% The following functions were originally part of the COBRA package. They
% are pasted here so that the step realizer can be distributed as a single
% m-file. Some of the code may be unreachable. 
% -------------------------------------------------------------------------
function Y = datahankel(d, m, rowdim, coldim)
%datahankel: Build a data Hankel matrix from a signal.
% 
%   Y = datahankel(D, m)
%
%       Construct a block-Hankel matrix of m block rows from the N-by-p
%       signal D, where is N is the number of samples in the signal and p
%       is the dimension of the signal.
% 
%   Y = datahankel(D, m, r, c)
%       
%       Construct a block-Hankel matrix of m block rows from the N-by-r*c
%       signal D, where N is the number of samples in the signal, r is the
%       row dimension of the signal, and c is the column dimension of the
%       signal. The first c columns of the i-th row of D will be the first
%       row of the i-th block of Y.
%       
%       If size(d, 3) > 1, then d is treated as an r x c x N signal. 
% 
%   The second version is intended for building data matrices with
%   matrix-valued signals, such as correlation functions. Because Matlab's
%   xcorr/xcov functions effectively vectorize signals, we must do some
%   reshaping to get them into block-matrix form.

% Uses custom functions
%   blkhankel

% (C) 2010 D. Miller, UCSD.
    error(nargchk(2, 4, nargin));

    if ~isscalar(m) || ~isreal(m)
        error('m must be real scalar.');
    end

    if nargin < 3
        if size(d, 3) > 1
            % If we have the r x c x N signal, then the dim arguments are
            % just taken from the signal dims.
            rowdim = size(d, 1);
            coldim = size(d, 2);
        else
            % Otherwise, we assume a p-dimensional signal.
            rowdim = size(d, 2);
            coldim = 1;
        end
    elseif nargin < 4
        error('Column dimension must be specified if row dimension is specified.');
    end

    % We compute the first block column of the data matrix and the last
    % block row, then let the blkhankel function do the rest of the work.
    if rowdim == 1 && coldim == 1 && size(d, 3) == 1
        % Dimension of the signal.
        dim = size(d, 2);

        % Put the signal into column-vector form.
        row = d';
        col = row(:);

        col = col(1:m*dim);
        row = row(:,m:end);
    else
        % Check if we've got a matrix signal, and reshape if so.
        if size(d, 3) > 1
            % In case the user supplied bad dimensions.
            if size(d, 1) ~= rowdim
                error('Invalid row dimension.');
            elseif size(d, 2) ~= coldim
                error('Invalid col dimension.');
            end
            d0 = d;
            d = zeros(size(d, 3), rowdim*coldim);
            for k = 1:rowdim
                for j = 1:coldim
                    d(:, (k-1)*coldim + j) = squeeze(d0(k, j, :));
                end
            end
        end

        N = size(d, 1);

        % Check that the dimensions make sense.
        if size(d, 2) ~= rowdim*coldim
            error('Number of columns of signal must = row*col.');
        end

        col = zeros(m*rowdim, coldim);
        for i = 1:m
            col(rowdim*(i-1)+1:rowdim*i, :) = reshape(d(i, :), coldim, rowdim)';
        end
        row = zeros(rowdim, coldim*(N-m+1));
        for i = 1:N-m+1
            row(:, coldim*(i-1)+1:coldim*i) = reshape(d(i+m-1, :), coldim, rowdim)';
        end
    end

    % Build a block-Hankel matrix from the data.
    Y = blkhankel(col, row);
end


function H = blkhankel(C, R)
%blkhankel: Construct a block-Hankel matrix.
% 
%   H = blkhankel(C, R) 
%
%       Construct a block-Hankel matrix with the first block colum C and
%       the last block row R. 
%
%       The total rows of C must be a multiple of the rows of R.

% (C) 2010 D. Miller

    error(nargchk(2, 2, nargin));

    cdim = size(C);
    rdim = size(R);

    blkrows = cdim(1)/rdim(1);
    blkcols = rdim(2)/cdim(2);

    % Error check dimensions of the input matrices.
    if mod(blkrows, 1) ~= 0
        error('Number of rows in C must be multiple of number of rows of R.');
    elseif mod(blkcols, 1) ~= 0
        error('Number of columns of R must be multiple of number of columns in R.');
    end

    coldim = cdim(2);
    rowdim = rdim(1);

    % Make sure the last block element of C and the first block element of
    % R are the same. No precedence rules as with the built-in hankel
    % function, just throw an error.
    if ~all(all(C(end-rowdim+1:end,:) == R(:,1:coldim)))
        error('SICLtools:blkhankel', ...
            'Last block-element of C must be first block-element of R.')
    end

    % Build the block-Hankel matrix.
    H = [C, zeros(cdim(1), rdim(2) - coldim)];
    colidx_prev = 1:coldim;
    for i = 2:blkcols
        colidx = (i-1)*coldim+1:i*coldim;
        H(:, colidx) = [H(rowdim+1:end, colidx_prev); R(:, colidx)];
        colidx_prev = colidx;
    end

end


function n = selectorder(s)
%selectorder: Select a system order from Hankel singular value estimates.
% 
%   n = selectorder(s)
% 
%       Return the selected order of the system with the singular values
%       contained in the vector s. The user selects a system order by 
%       clicking on the corresponding value in a plot of singular values.
% 
%       If the user exits without selecting an order, then n = -1.

% (C) 2010 D. Miller, UCSD.
    
    if nargin ~= 1 || ~isvector(s) || ~isreal(s)
        error('Only one real vector argument permitted.');
    end

    % Create a figure and move the axis up from the default slightly.
    hFigure = figure;
    figPos = get(hFigure, 'Position');
    figColor = get(hFigure, 'Color');
    set(hFigure, 'Position', [figPos(1), figPos(2), 560, 420]);
    hAxis = axes('Parent', hFigure, 'Units', 'Pixels', ...
        'Position', [75 55 433 335]);

    % Plot the singular values.
    semilogy(hAxis, 1:length(s), s, 'ob', 'MarkerFaceColor', 'blue', 'MarkerSize', 10);
    title('Choose system order...')
    ylabel('Singular Value')
    axis([0.9, length(s)+0.1, min(s), max(s)]);

    % Add a text box on the bottom to show what order the user selected.
    basePos = [75, 5];
    orderWidth = 95;
    uicontrol('Style', 'text', 'String', 'Selected Order: ', ...
        'Position', [basePos(1), basePos(2), orderWidth, 20], ...
        'BackgroundColor', figColor);
    hOrderSelection = uicontrol('Style', 'text', 'String', '0', ...
        'Position', [basePos(1)+orderWidth, basePos(2), 30, 20], ...
        'BackgroundColor', figColor);


    % Create a local-scope callback function to update the system order.
        n = -1;
        function update_n(k)
            n = k;
        end

    % Add a callback to the figure window to detect when the user clicks on
    % the plot.
    axis_callback = @(h, x) plot_callback(h, x, s, ...
        hAxis, hOrderSelection, @update_n);
    set(hFigure, 'WindowButtonDownFcn', axis_callback);


    % Not the most elegant way of doing this, but Matlab's OO stuff doesn't
    % seem to apply to figures... We pause for 0.1 sec until exitFunc becomes
    % true. This happens when the button is pressed, which runs the callback
    % exit_callback.
        exitFunc = false;
        function exit_callback(~, ~)
            exitFunc = true;
            delete(hFigure);
        end

    uicontrol('Style', 'pushbutton', 'String', 'OK', ...
        'Position', [basePos(1)+orderWidth+280, basePos(2)+5, 60, 20], ...
        'Callback', @exit_callback);

    % We also want to exit if the user closes the window.
    set(hFigure, 'CloseRequestFcn', @exit_callback);

    pause on;
    while ~exitFunc
        pause(0.1);
    end
end

function plot_callback(~, ~, s, hAxis, hOrderSelection, update_n)
    point = get(hAxis, 'CurrentPoint');

    % Ignore clicks outside the axis window.
    if point(1, 2) > max(s) || point(1, 2) < min(s) || ...
       point(1, 1) < 0.9 || point(1, 1) > length(s)+0.1
        return;
    else
        n = round(point(1,1));
    end

    % Redraw all dots, then redraw the selected dot as red.
    line(1:length(s), s, ...
        'LineStyle', 'none', ...
        'Color', 'blue', ...
        'Marker', 'o', ...
        'MarkerFaceColor', 'blue', ...
        'MarkerSize', 10);
    line(n, s(n),  ...
        'LineStyle', 'none', ...
        'Color', 'red', ...
        'Marker', 'o', ...
        'MarkerFaceColor', 'red', ...
        'MarkerSize', 10);
    
    % Update order selection text box.
    set(hOrderSelection, 'String', num2str(n));
    
    update_n(n);
end
