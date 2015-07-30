classdef Dynamics < dscomponents.ACompEvalCoreFun
    % This class implements the nonlinear continuum mechanics as described
    % in @cite Heidlauf2013 .
    
    properties
        % Cross-fibre markert part
        b1cf = 5.316372204148964; % [MPa] = [N/mm²]
        d1cf = 0.014991843974911; % [-]
        
        % The activation of the muscle at time t
        %
        % @type function_handle @default @(t)0
        alpha = @(t)0; % [-]
        
        %% Unassembled stuff
        ComputeUnassembled = false;
        % Sigma assembly matrix
        Sigma;
        % The indices of any dirichlet value in the unassembled vector duvw
        idx_uv_bc_glob_unass;
        NumTotalDofs_unass;
        idx_vp_dof_unass_elems;
        
        ForceLengthFun;
        ForceLengthFunDeriv;
    end
    
    properties(SetAccess=private)
        APExp;
        AnisoPassiveMuscle;
        AnisoPassiveMuscleDeriv;
        AnisoPassiveTendon;
        AnisoPassiveTendonDeriv;
    end
    
    properties(SetAccess=protected)
        % Helper variable for fullmuscle.model
        lambda_dot_pos;
        lambda_dot;
        
        nfibres;
    end
    
    properties(Transient, SetAccess=private)
        % Prepared arguments for APExpansion
        %         muprep;
        LastBCResiduals;
        
        nfevals;
        nJevals;
    end
    
    properties(Transient)
        % Helper value for QuickReleaseTests (or others) that use a
        % function handle with certain alpha ramp time for different
        % simulations. Used in getSimCacheExtra to uniquely identify a
        % simulation in the cache
        RampTime;
        
        % Cached values of the current constraint evaluation and it's
        % jacobian
        curGC;
        curJGC;
    end
    
    properties(Transient, Access=private)
        % Cached quantity from this.fsys.UseDirectMassInversion for
        % faster evaluation of dynamics.
        usemassinv;
        
        % Cached value for cross fibre computations (speed)
        crossfibres = false;
    end
    
    properties(Access=private)
        % Reference to the full system
        fsys;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            
            %% Load AP expansion
            %             d = fileparts(which('muscle.Dynamics'));
            %             s = load(fullfile(d,'AP'));
            %             s.kexp.Ma = s.kexp.Ma(1,:);
            %             this.APExp = s.kexp;
        end
        
        function configUpdated(this)
            sys = this.fsys;
            mc = sys.Model.Config;
            if ~isempty(mc.FibreTypeWeights)
                this.nfibres = size(mc.FibreTypeWeights,2);
            end
            if ~isempty(mc)
                this.xDim = sys.NumTotalDofs;
                this.fDim = sys.NumDerivativeDofs;
                this.JSparsityPattern = this.computeSparsityPattern;
                
                %% Sigma assembly matrix
                this.precomputeUnassembledData;
            end
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACompEvalCoreFun(this, mu);
            
            % Reset counters
            this.nfevals = 0;
            this.nJevals = 0;
            
            sys = this.fsys;
            mc = sys.Model.Config;
            if ~isempty(mc.Pool)
                mc.Pool.prepareSimulation(sys.Model.T,mu(4));
            end
            % Returns an all zero function if mu(2) is less or equal to zero!
            this.alpha = mc.getAlphaRamp(mu(2));
            
            % Prepare force-length fun
            mc.setForceLengthFun(this);
            
            % Muscle anisotropic passive law
            mlfg = tools.MarkertLawOriginal(mu(5),mu(6));
            [this.AnisoPassiveMuscle, this.AnisoPassiveMuscleDeriv] = mlfg.getFunction;
            
            % Tendon anisotropic passive law
            mlfg = tools.CubicToLinear(mu(7),mu(8));
%             mlfg = tools.MarkertLaw(1.3895e+07,11.1429,1.637893706954065e+05);
            [this.AnisoPassiveTendon, this.AnisoPassiveTendonDeriv] = mlfg.getFunction;

            % Get the law function handles that also take b,d as arguments.
            % Needed due to possibly inhomogeneous material.
            %[~,~,this.MarkertLawFun,this.MarkertLawFunDeriv] = mlfg.getFunction;
            
            % Cache stuff
            this.usemassinv = sys.UseDirectMassInversion;
            this.crossfibres = sys.HasFibres && sys.UseCrossFibreStiffness;
        end
        
        function fx = evaluateCoreFun(this, x, t)
            error('Do not call directly; have custom evaluate method.');
        end
        
        function fx = evaluateComponentSet(this, nr, x, t)
            % Computes the full or reduced component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            %
            % See also: PointSet
            fx = this.evaluate(x,t);
            fx = this.V(this.PointSets{nr},:)*fx;
        end
        
        function fx = evaluateComponentSetMulti(this, nr, x, t, mu)
            % Computes the full component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space locations `\vx` @type matrix<double>
            % t: The corresponding times `t` for the state `\vx` @type
            % rowvec<double>
            % mu: The corresponding parameters `\mu`. Can be a single
            % parameter or as many as the size of x @type matrix<double>
            %
            % See also: PointSet
            fx = this.evaluateMulti(x,t,mu);
            fx = fx(this.PointSets{nr},:);
        end
        
        function res = test_Jacobian(this, y, t, mu)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            %
            % Currently the tests using viscosity are commented out as we
            % assume linear damping, which is extracted as extra `A(t,\mu)`
            % part in the models' system
            
            if nargin < 4
                mu = this.System.Model.DefaultMu;
                if nargin < 3
                    t = 1000;
                    if nargin < 2
                        y = this.System.x0.evaluate(mu);
                    end
                end
            end
            
            % Use nonzero t to have an effect
            res = test_Jacobian@dscomponents.ACoreFun(this, y, t, mu);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = muscle.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far everything is done in the
            % constructor.
        end
        
        function setSystem(this, sys)
            setSystem@dscomponents.ACoreFun(this, sys);
            if isa(sys,'models.ReducedSecondOrderSystem')
                this.fsys = sys.Model.FullModel.System;
            else
                this.fsys = sys;
            end
        end
        
        function copy = project(this, V, W)
             copy = project@dscomponents.ACompEvalCoreFun(this, V, W);
        end
        
        function set.ComputeUnassembled(this, value)
            if value && ~isempty(this.NumTotalDofs_unass)%#ok
                this.fDim = this.NumTotalDofs_unass;%#ok
            elseif ~value && ~isempty(this.System)
                this.fDim = this.System.NumTotalDofs;
            end
            this.ComputeUnassembled = value;
        end
        
        % Declare public
        [SPK, SPg, SPalpha, SPLamDot] = computeSparsityPattern(this);
    end
    
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, ~, ~, x, ~)
            % This is the template method that actually evaluates the components at given points
            % and values.
            %
            % @attention This method must be able to handle vector-arguments
            % for `\vx,t,\vmu`!
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1)+1{:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: A matrix `\vX` with the state space locations `\vx_i` in
            % its columns. In rows end(i-1)+1:end(i) it contains the states
            % relevant to evaluate the i-th component of �\vf�.
            % States occur multiply in �\vx� if different components
            % of �\vf� depend on these states.
            % @type matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, dfxsel)%#ok
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            % mu: The corresponding parameter `\mu` for the state `\vx` @type colvec<double>
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            
        end
    end
    
    methods(Access=private)
        function precomputeUnassembledData(this)
            sys = this.System;
            mc = sys.Model.Config;
            geo = mc.PosFE.Geometry;
            num_u_glob = 3*geo.NumNodes;
            outsize = num_u_glob;
            
            % Position part: not assembly as u' = v without FEM
            
            % Velocity part: x,y,z velocities
            [i, ~] = find(mc.PosFE.Sigma);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            
            % Pressure part
            pgeo = mc.PressFE.Geometry;
            [i, ~] = find(mc.PressFE.Sigma);
            I = [I(:); outsize+i];
            outsize = outsize + pgeo.NumNodes;
            
            n = numel(I);
            S = sparse(I,1:n,ones(n,1),outsize,n);
            
            % Take out nodes with dirichlet BC on output side
            S([sys.idx_u_bc_glob; sys.idx_v_bc_glob-num_u_glob],:) = [];
            % Find corresponding unassembled dofs that would be ignored
            % (due to dirichlet velocity values, pressure dirichlet not
            % implemented)
            bc_unass = find(sum(S,1) == 0);
            % Remove them, too. The unassembled evaluation also removes the
            % corresponding entries of the unassembled vector.
            S(:,bc_unass) = [];
            this.idx_uv_bc_glob_unass = [sys.idx_u_bc_glob' num_u_glob + bc_unass];
            this.Sigma = S;
            
            this.NumTotalDofs_unass = sys.NumStateDofs + size(S,2);
            
            % Create boolean array that determines which unassembled dofs
            % belong to which element
            % dv dofs
            hlp = repmat(1:geo.NumElements,3*geo.DofsPerElement,1);
            pgeo = mc.PressFE.Geometry;
            % dp dofs
            hlp2 = repmat(1:geo.NumElements,pgeo.DofsPerElement,1);
            hlp = [hlp(:); hlp2(:)];
            hlp(bc_unass) = [];
            ass = false(geo.NumElements,length(hlp));
            for k = 1:geo.NumElements
                ass(k,:) = hlp == k;
            end
            this.idx_vp_dof_unass_elems = ass;
        end
    end
    
    methods(Static)
        function res = test_UnassembledEvaluation
            m = muscle.Model(Cube12);
            mu = m.DefaultMu;
            [t,~,~,x] = m.simulate(mu);
            s = m.System;
            f = s.f;
            fx = f.evaluate(x(:,1),t(1));
            f.ComputeUnassembled = true;
            fxu = f.evaluate(x(:,1),t(1));
            fx_ass = zeros(size(fx,1),1);
            % du part (no assembly)
            fx_ass(1:s.NumStateDofs) = fxu(1:s.NumStateDofs);
            % dvw part (with assembly)
            fx_ass(s.NumStateDofs + (1:s.NumDerivativeDofs+s.NumAlgebraicDofs)) = s.f.Sigma * fxu(s.NumStateDofs+1:end);
            res = Norm.L2(fx - fx_ass) < eps;
            
            % Multi-eval
            f.ComputeUnassembled = false;
            fx = f.evaluateMulti(x,t,mu);
            f.ComputeUnassembled = true;
            fxu = f.evaluateMulti(x,t,mu);
            % du part (no assembly)
            fx_ass = zeros(size(fx,1),length(t));
            fx_ass(1:s.NumStateDofs,:) = fxu(1:s.NumStateDofs,:);
            % dvw part (with assembly)
            fx_ass(s.NumStateDofs + (1:s.NumDerivativeDofs+s.NumAlgebraicDofs),:) = s.f.Sigma * fxu(s.NumStateDofs+1:end,:);
            res = res & sum(Norm.L2(fx - fx_ass)) < eps;
        end
    end
end
