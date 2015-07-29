classdef System < models.BaseSecondOrderSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% The system u'' + K(u) is transformed into the first order system v' +
% K(u), u' = v
%
% @author Daniel Wirtz @date 2014-01-20

    properties
       % The global index of node x,y,z positions within uvw.
       %
       % The velocities are indentically indexed but 3*NumNodes later.
       idx_u_glob_elems;
       
       % The global index of node pressures within uvw.
       idx_p_glob_elems;
    end
    
    properties(SetAccess=private)
       % This contains the indices of the nodes tracking pressure (linear /
       % 8-corner elements) within the quadratic global node numbering
       %
       % Used for plotting only so far.
       idx_p_to_u_nodes;
       
       % The overall values of dirichlet conditions
       %
       % To specify in [mm] for position and [mm/ms] for velocity
       % 
       % Collected from val_u_bc, val_v_bc
       val_uv_bc_glob;
       
       % The positions of the dirichlet values within the global node/xyz
       % state space vector
       %
       % Collected from idx_u_bc_glob, idx_v_bc_glob
       idx_uv_bc_glob;
       
       %% Position dirichlet fields
       % Boundary conditions: 3 times numnodes logical matrix telling
       % which degree of freedom of which (position) node is fixed
       %     n1   n2   n3 ...
       % x | 1    0    0       |
       % y | 1    1    0 ...   |
       % z | 1    0    0       |
       bool_u_bc_nodes;
       idx_u_bc_glob;
       idx_u_bc_local;
       val_u_bc; % [mm]
       
       %% Velocity dirichlet fields
       % This set comprises the explicitly/directly set velocity boundary
       % conditions
       bool_expl_v_bc_nodes;
       idx_expl_v_bc_glob;
       idx_expl_v_bc_local;
       val_expl_v_bc; % [mm/ms]
       
       % This set comprises both the explicit velocity conditions and the
       % implicit zero velocity conditions that apply for each position
       % dirichlet condition
       idx_v_bc_glob;
       idx_v_bc_local;
       val_v_bc; % [mm/ms]
       
       % Positions of velocity-bc affected DoFs in the u dofs
       idx_v_bc_u_dof;
              
       %% Neumann fields
       bc_neum_forces_nodeidx; % [N]
       bc_neum_forces_val;
       FacesWithForce;
       
       %% Convenience variables
       % A helper array containing the indices of the actual dofs in the
       % global indexing.
       %
       % Used to join dofs with dirichlet values
       idx_uv_dof_glob;
       
       % The indices of u,v,p components in the effective dof vector
       idx_u_dof_glob;
       idx_v_dof_glob;
       num_v_bc;
       idx_p_dof_glob;
       
       % Total number of discrete field points including dirichlet values
       num_uvp_glob;
       
       Minv;
       
       %% Fibre stuff
       HasFibres = false;
       HasFibreTypes = false;
       HasMotoPool = false;
       % Property only useful when fullmodel.System is used.. comes due to
       % use of inheritance. Better solutions could be thought of, but this
       % is to get it going!
       HasForceArgument = false;
       a0;
       a0oa0;
       dNa0;
       a0Base;
       
       %% Tendon stuff
       HasTendons = false;
       % The function used to map between muscle (mi) and tendon (ma) parameter over
       % [0,1] (v)
       MuscleTendonParamMapFun = @(v,mi,ma)10.^(log10(mi) + v.*(log10(ma)-log10(mi)));
       MuscleTendonRatioGP = [];
       MuscleTendonRatioNodes = [];
       
       % Fields to contain c10/c01 values for mooney-rivlin law
       % (muscle+tendon)
       MuscleTendonParamc10 = [];
       MuscleTendonParamc01 = [];
       % Stress free initial condition constant for mooney-rivlin law
       % This was "classically" done by setting the pressure to this
       % constant, however, if inhomogeneous material constants are
       % present, the linear pressure model cannot resolve those values and
       % hence any initial condition would be inconsistent.
       % Now, the pressure is initialized as zero and the
       % MooneyRivlinICConst takes care of consistent initial conditions.
       MooneyRivlinICConst = [];
       
       
       %% Cross-fibre stiffness stuff
       % normals to a0
       a0oa0n1;
       a0oa0n2;
    end
   
    properties(Dependent)
       % Flag to invert the velocity mass matrix before simulations.
       %
       % If this is set to true, the mass matrix for the velocity part will
       % be assembled, the boundary condition rows removed and the inverse
       % of the remaining matrix will be pre-computed.
       % This inverse will be pre-multiplied to the velocity-dofs inside
       % the muscle.Dynamics evaluate and getStateJacobian functions.
       %
       % This in general results in higher simulation speed but it remains
       % open to see how well reduced modeling will work with that scheme.
       %
       % @type logical @default false
       UseDirectMassInversion;
       
       % Flag to indicate if this system should use stiffness terms
       % (Markert) for cross-fibre directions.
       %
       % @type logaical @default false
       UseCrossFibreStiffness;
    end
    
    properties(Access=private)
        fUseDirectMassInversion = false;
        fUseCrossFibreStiffness = false;
        fD;
    end
    
    properties(SetAccess=protected)
        Plotter;
    end
    
    properties(Transient, Access=private)
        % Cached velocity bc time-dependent function
        velo_bc_fun = [];
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseSecondOrderSystem(model);
            
            % The muscle viscosity
            % @unit [g / (mm * ms)] = [kP] (kiloPoiseulle)
            this.addParam('viscosity',.1);
            
            % The amount of milliseconds over which to activate the model
            % Set this parameter to zero to deactivate (any maybe use a
            % custom activation ramp)
            %
            % @unit [ms]
            this.addParam('alpha_ramp_time',0);
            
            % The pressure on the faces (neumann conditions)
            %
            % @unit [MPa] = [N / mm^2]
            this.addParam('neumann_pressure',0);
            
            % For some variants, we have the mean input current for the
            % motoneuron pool (generating activation)
            this.addParam('mean input current',0);
            
            % anisotropic passive stiffness for muscle material
            % markert law b1
            %
            % @unit [MPa]
            % #5
            this.addParam('muscle passive b1',2.756e-5);
            
            % anisotropic passive stiffness for muscle material
            % markert law d1 [-]
            % #6
            this.addParam('muscle passive d1',43.373);
            
            % anisotropic passive stiffness for tendon material
            % markert law
            %
            % @unit [MPa]
            %
            % -or-
            % ideal length of sarcomere w.r.t force production [-]
            %
            % fit with tools.MarkertLaw: b1 = 1.3895e+04
            % fit with tools.QuadToLinear: lam0 = 1.027
            % fit with tools.CubicToLinear: lam0 =  1.0175 [26.2.15]
            % fit with tools.CubicToLinear: lam0 =  1.0207 [29.4.15, after scaling fix]
            % #7
            this.addParam('tendon passive 1 [b1/lam0]', 1.0207); % [MPa/-]
            
            % anisotropic passive stiffness for tendon material
            % markert law
            %
            % fit with tools.MarkertLaw: d1 = 11.1429 [-]
            % max modulus  1.637893706954065e+05
            % fit with tools.QuadToLinear: M = 1.637893706954065e+05
            % fit with tools.CubicToLinear: M =  1.637893706954065e+05 [26.2.15]
            % fit with tools.CubicToLinear: M =  1.637893706954065e+05 [29.4.15, after scaling fix]
            % #8
            this.addParam('tendon passive 2 [d1/M]', 1.637893706954065e+05);
            
            % isotropic muscle material
            % mooney-rivlin law c10
            %
            % from sprenger thesis: 35.6 kPa
            %
            % @unit [MPa]
            %
            % #9
            this.addParam('muscle mooney-rivlin c10',35.6e-3);
            
            % isotropic muscle material
            % mooney-rivlin law c01
            %
            % from sprenger thesis: 3.86 kPa
            %
            % @unit [MPa]
            %
            % #10
            this.addParam('muscle mooney-rivlin c01',3.86e-3);
            
            % isotropic tendon material
            % mooney-rivlin law c10
            %
            % @unit [MPa]
            %
            % #11
            this.addParam('tendon mooney-rivlin c10',2.310e3);
            
            % isotropic tendon material
            % mooney-rivlin law c01
            %
            % @unit [MPa]
            %
            % #12
            this.addParam('tendon mooney-rivlin c01',1.15e-3);
            
            % muscle fibre maximal force
            %
            % @unit [MPa]
            %
            % #13
            this.addParam('P_max',.3);
            
            % parameter p1 for force-length curve customization.
            %
            % For linear curve from Gordon66 we set lambda_0 with this,
            % i.e. the resting sarcomere length. According to literature,
            % somewhere between 2 and 2.2 micrometer.
            % #14
            this.addParam('force-length p1 (lam_0/width/...)',2.05);
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
        end
        
        function rsys = buildReducedSystem(this, rmodel)
            % Overrides the default buildReducedModel method in order to
            % temporarily set the A component used in projection.
            %
            % (So far the extra efficiency of "no A" for mu(1) == 0 is
            % neglected in any reduced model)
            this.D = this.fD;
            rsys = buildReducedSystem@models.BaseSecondOrderSystem(this, rmodel);
            this.D = [];
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            if ~isempty(mc)
                tq = mc.PosFE;
                geo_pos = tq.Geometry;
                tl = mc.PressFE;
                geo_p = tl.Geometry;
                
                % Find the indices of the pressure nodes in the displacement
                % nodes geometry (used for plotting)
                this.idx_p_to_u_nodes = geo_pos.getCommonNodesWith(geo_p);

                % Call subroutine for boundary condition index crunching
                % This needs to be done before we can compute the effective
                % Dofs of the system.
                this.computeDirichletBC;
                
                %% Get dimensions
                this.updateDimensions(mc);
                
                % Construct B matrix
                this.B = this.assembleB;
                
                % Set input function
                this.Inputs = mc.getInputs;

                %% Tendon stuff
                % Detect of tendons are present
                this.initMuscleTendonRatios;
                
                % Init fibre directions and precomputable values
                this.inita0;
                
                this.HasMotoPool = this.HasFibres && ~isempty(mc.Pool);
                this.HasFibreTypes = this.HasFibres && ~isempty(mc.FibreTypeWeights);
                this.HasForceArgument = this.HasFibreTypes && isa(this,'fullmuscle.System');

                % Construct global indices in uw from element nodes. Each dof in
                % an element is used three times for x,y,z displacement. The
                % "elems" matrix contains the overall DOF numbers of each
                % element in the order of the nodes (along row) in the master
                % element.
                ne = geo_pos.NumElements;
                globalelementdofs = zeros(3,geo_pos.DofsPerElement,ne,'int32');
                for m = 1:ne
                    % First index of element dof in global array
                    hlp = (geo_pos.Elements(m,:)-1)*3+1;
                    % Use first, second and third as positions.
                    globalelementdofs(:,:,m) = [hlp; hlp+1; hlp+2];
                end
                this.idx_u_glob_elems = globalelementdofs;

                % The same for the pressure
                globalpressuredofs = zeros(geo_p.DofsPerElement,geo_p.NumElements,'int32');
                off = 2 * (geo_pos.NumNodes * 3);
                for m = 1:geo_p.NumElements
                    globalpressuredofs(:,m) = off + geo_p.Elements(m,:);
                end
                this.idx_p_glob_elems = globalpressuredofs;

                %% Compile Mass Matrix
                this.M = dscomponents.ConstMassMatrix(this.assembleMassMatrix);

                %% Compile Damping Matrix
                this.fD = this.assembleDampingMatrix;

                %% Tell f we have a new config
                this.f.configUpdated;
                
                %% Initial value
                this.x0 = dscomponents.ConstInitialValue(this.assembleX0);
                
                %% Algebraic constraints function
                this.g = muscle.ConstraintsFun(this);
                
                this.updateSparsityPattern;
            end
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Check for viscosity setting
            %
            % See also: muslce.System.MooneyRivlinICConst
            
            this.D = [];
            if mu(1) > 0
                this.D = this.fD;
            end
            % Update muscle/tendon parameters on all gauss points
            % We have mu-dependent mooney-rivlin and markert laws, but they
            % are constant over one simulation. So precompute here.
            this.updateTendonMuscleParamsGP(mu);
            
            % Update the MooneyRivlinICConst to have stress-free IC
            % This depends on the current muscle/tendon parameters for
            % mooney-rivlin (updated one step before)
            this.MooneyRivlinICConst = -2*this.MuscleTendonParamc10...
                -4*this.MuscleTendonParamc01;
            
            % Get velocity dirichlet conditions time-function
            mc = this.Model.Config;
            if ~isempty(mc.VelocityBCTimeFun)
                this.velo_bc_fun = mc.VelocityBCTimeFun.getFunction;
            end
            
            prepareSimulation@models.BaseSecondOrderSystem(this, mu, inputidx);
        end
        
        function pm = plotDiff(this, t, uvw1, uvw2, fac, varargin)
            if nargin < 5
                fac = 5;
            end
            x0 = this.x0.evaluate([]);
            diff = repmat(x0,1,length(t)) + (uvw1-uvw2)*fac;
            pm = this.plot(t,diff,varargin{:});
        end
        
        function uvwall = includeDirichletValues(this, t, uvw)
            %% Re-add the dirichlet nodes
            % Efficient for single vector
            if size(uvw,2) == 1
                uvwall = zeros(this.num_uvp_glob, 1);
                uvwall(this.idx_uv_dof_glob) = uvw;
                uvwall(this.idx_uv_bc_glob) = this.val_uv_bc_glob;
                if ~isempty(this.velo_bc_fun)
                    uvwall(this.idx_v_bc_glob) = uvwall(this.idx_v_bc_glob)*this.velo_bc_fun(t);
                end
            else
                uvwall = zeros(this.num_uvp_glob, size(uvw,2));
                uvwall(this.idx_uv_dof_glob,:) = uvw;
                uvwall(this.idx_uv_bc_glob,:) = repmat(this.val_uv_bc_glob,1,size(uvw,2));
                if ~isempty(this.velo_bc_fun)
                    uvwall(this.idx_v_bc_glob,:) = bsxfun(@times, uvwall(this.idx_v_bc_glob,:), this.velo_bc_fun(t));
                end
            end
        end
    end
    
    methods
        function set.UseDirectMassInversion(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('UseDirectMassInversion must be true or false');
            end
            if this.fUseDirectMassInversion ~= value
                this.fUseDirectMassInversion = value;
                this.configUpdated;
            end
        end
        
        function value = get.UseDirectMassInversion(this)
            value = this.fUseDirectMassInversion;
        end
        
        function set.UseCrossFibreStiffness(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('UseCrossFibreStiffness must be true or false');
            end
            if this.fUseCrossFibreStiffness ~= value
                this.fUseCrossFibreStiffness = value;
                if ~isempty(this.Model.Config)
                    this.inita0;
                end
            end
        end
        
        function value = get.UseCrossFibreStiffness(this)
            value = this.fUseCrossFibreStiffness;
        end
    end
    
    methods(Access=protected)
        
        function updateDimensions(this, mc)
            tl = mc.PressFE;
            geo_p = tl.Geometry;
            this.NumAlgebraicDofs = geo_p.NumNodes;
            
            tq = mc.PosFE;
            geo_uv = tq.Geometry;
            this.NumStateDofs = geo_uv.NumNodes * 3 - length(this.idx_u_bc_glob);
            
            this.num_uvp_glob = geo_uv.NumNodes * 6 + this.NumAlgebraicDofs;
            % Compile the inverse addressing for inserting dofs into the
            % full vector
            idx = int32(1:this.num_uvp_glob);
            idx(this.idx_uv_bc_glob) = [];
            this.idx_uv_dof_glob = idx;
            
            % Have same number of xdot dofs than x dofs except the case
            % that explicit velocity dirichlet conditions are provided
            this.NumDerivativeDofs = this.NumStateDofs - length(this.idx_expl_v_bc_local);
            
            updateDimensions@models.BaseSecondOrderSystem(this);
        end
        
        function x0 = assembleX0(this)
            % Constant initial values as current node positions
            mc = this.Model.Config;
            tq = mc.PosFE;
            geo = tq.Geometry;
            
            % Fill in the reference configuration positions as initial
            % conditions
            x0 = geo.Nodes(:);
            % Remove those fixed via dirichlet conditions
            x0(this.idx_u_bc_local) = []; %local=global here as is first set
            
            % Give the model config a chance to mess with x0
            x0 = mc.getX0(x0);
        end
        
        function Baff = assembleB(this)
            % Collect neumann forces
            [B, this.bc_neum_forces_nodeidx] = this.getSpatialExternalForces;
            % Only set up forces if present
            Baff = [];
            if ~isempty(this.bc_neum_forces_nodeidx)
                this.bc_neum_forces_val = B(this.bc_neum_forces_nodeidx);
                % Remove dirichlet DoFs
                B(this.idx_v_bc_local) = [];
                Baff = dscomponents.AffLinInputConv;
                Baff.TimeDependent = false;
                Baff.addMatrix('mu(3)',B);
            end
        end
        
        function MM = assembleMassMatrix(this)
            %% Compile Mass Matrix
            mc = this.Model.Config;
            fe_pos = mc.PosFE;
            g = fe_pos.Geometry;
            
            % Augment mass matrix for all 3 displacement directions
            nd = g.NumNodes;
            [i, j, s] = find(fe_pos.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            MM = this.Model.MuscleDensity*sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            
            % Strip out the entries of dirichlet nodes
            MM(this.idx_v_bc_local,:) = [];
            MM(:,this.idx_v_bc_local) = [];
            
            % See description of property
            if this.UseDirectMassInversion
                this.Minv = inv(MM);
                MM = 1;
            end
        end
        
        function Daff = assembleDampingMatrix(this)
            %% Compile Damping/Viscosity Matrix
            mc = this.Model.Config;
            fe_pos = mc.PosFE;
            g = fe_pos.Geometry;
            
            % Augment mass matrix for all 3 displacement directions
            nd = g.NumNodes;
            [i, j, s] = find(fe_pos.D);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            D = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            
            % Strip out the entries of dirichlet nodes
            D(this.idx_v_bc_local,:) = [];
            D(:,this.idx_v_bc_local) = [];
            
            Daff = dscomponents.AffLinCoreFun(this);
            Daff.addMatrix('mu(1)',-D);
        end
        
        function computeDirichletBC(this)
            mc = this.Model.Config;
            [pos_dir, velo_dir, velo_dir_val] = mc.getBC;
            
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            fe_press = mc.PressFE;
            pgeo = fe_press.Geometry;
            
            % Total number of position entries in global vector
            num_u_glob = geo.NumNodes * 3;
            
            %% Displacement / State Dirichlet
            this.bool_u_bc_nodes = pos_dir;
            % Set values to node positions
            this.val_u_bc = geo.Nodes(pos_dir);
            this.idx_u_bc_local = int32(find(pos_dir(:)));
            % Local = global as u is the first set of variables in the
            % global vector
            this.idx_u_bc_glob = this.idx_u_bc_local;
            
            %% Velocity / Derivative Dirichlet
            % Add any user-defines values
            this.bool_expl_v_bc_nodes = velo_dir;
            this.idx_expl_v_bc_local = int32(find(velo_dir(:)));
            this.val_expl_v_bc = velo_dir_val(velo_dir);
            this.idx_expl_v_bc_glob = this.idx_expl_v_bc_local + num_u_glob;
            
            % Compute position of explicit derivative dirichlet conditions
            % inside the state dofs
            pos = false(num_u_glob,1);
            % Mark positions of explicit conditions
            pos(this.idx_expl_v_bc_local) = true;
            % Remove implicit velocity dirichlet conditions from position
            % dirichlet conditions
            pos(this.idx_u_bc_local) = [];
            this.DerivativeDirichletPosInStateDofs = pos;
            
            % Include the zero velocity conditions that apply for each
            % position dirichlet condition
            % (cannot conflict with position dirichlet conditions, this is
            % checked in AModelConfig.getBC)
            this.idx_v_bc_local = [this.idx_expl_v_bc_local; this.idx_u_bc_local];
            % Here we add zero velocities for each point with fixed position, too.
            this.val_v_bc = [this.val_expl_v_bc; zeros(size(this.val_u_bc))];
            this.num_v_bc = length(this.val_v_bc);
            this.idx_v_bc_glob = this.idx_v_bc_local + num_u_glob;
            
            %% Convenience index sets
            % Compile the global dirichlet values index and value vectors.
            this.idx_uv_bc_glob = [this.idx_u_bc_glob; this.idx_v_bc_glob];
            this.val_uv_bc_glob = [this.val_u_bc; this.val_v_bc];
            
            % Compute dof positions in global state space vector
%             
%             pos = false(1,total);
%             pos(1:num_u_glob) = true;
%             pos(this.idx_uv_bc_glob) = [];
%             this.idx_u_dof_glob = int32(find(pos));
            
%             pos = false(1,total);
%             pos(num_u_glob+1:num_u_glob*2) = true;
%             pos(this.idx_uv_bc_glob) = [];
%             this.idx_v_dof_glob = int32(find(pos));
            
%             pos = false(1,total);
%             pos(this.idx_v_bc_glob-num_u_glob) = true;
%             pos(this.idx_u_bc_glob) = [];
%             this.idx_v_bc_u_dof = int32(find(pos));
            
            % no possible dirichlet values for p. so have all
%             this.idx_p_dof_glob = geo.NumNodes*6-length(this.idx_uv_bc_glob) + (1:pgeo.NumNodes);
%             this.NumAlgebraicDofs = length(this.idx_p_dof_glob);
            % Automatically mark the pressure DoFs as algebraic condition.
            % Important for reduction.
%             this.AlgebraicConditionDoF = this.idx_p_dof_glob;
            
            
        end
        
        function val = getDerivativeDirichletValues(this, t)
            % Computes the derivative dirichlet values dependent on the
            % current time.
            %
            % See also: ODEFun DerivativeDirichletPosInStateDofs

            % Check if velocity bc's should be applied in time-dependent manner
            if ~isempty(this.velo_bc_fun)
                val = this.velo_bc_fun(t)*this.val_expl_v_bc;
            else
                val = this.val_expl_v_bc;
            end
        end
        
        function [force, nodeidx] = getSpatialExternalForces(this)
            mc = this.Model.Config;
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            ngp = fe_displ.GaussPointsPerElemFace;
            force = zeros(geo.NumNodes * 3,1);
            faceswithforce = false(1,geo.NumFaces);
            
            globalcoord = strcmp(mc.NeumannCoordinateSystem,'global');
            for fn = 1:geo.NumFaces
                elemidx = geo.Faces(1,fn);
                faceidx = geo.Faces(2,fn);
                masterfacenodeidx = geo.MasterFaces(faceidx,:);
                % So far: Constant pressure on all gauss points!
                P = mc.getBoundaryPressure(elemidx, faceidx);
                if ~isempty(P)
                    faceswithforce(fn) = true;
                    integrand = zeros(3,geo.NodesPerFace);
                    if globalcoord
                        N = geo.FaceNormals(:,faceidx);
                    end
                    for gi = 1:ngp
                        if ~globalcoord
                            N = fe_displ.NormalsOnFaceGP(:,gi,fn);
                        end
                        PN = (P * N) * fe_displ.Ngpface(:,gi,fn)';
                        integrand = integrand + fe_displ.FaceGaussWeights(gi)*PN*fe_displ.face_detjac(fn,gi);
                    end
                    facenodeidx = geo.Elements(elemidx,masterfacenodeidx);
                    facenodeidx = (facenodeidx-1)*3+1;
                    facenodeidx = [facenodeidx; facenodeidx+1; facenodeidx+2];%#ok
                    force(facenodeidx(:)) = force(facenodeidx(:)) + integrand(:);
                end
            end
            % Augment to u,v,w vector
            nodeidx = find(abs(force) > 1e-13);
            this.FacesWithForce = faceswithforce;
        end
        
        function inita0(this)
            mc = this.Model.Config;
            fe = mc.PosFE;
            geo = fe.Geometry;
            
            anull = mc.geta0;
            ismastercoord = strcmp(mc.a0CoordinateSystem,'master');
            this.HasFibres = false;
            if any(anull(:))
                this.HasFibres = true;
            
                this.a0 = anull;

                % Precomputations
                dNgp = fe.gradN(fe.GaussPoints);
                Ngp = fe.N(fe.GaussPoints);
                
                anulldyadanull = zeros(3,3,fe.GaussPointsPerElem*geo.NumElements);
                dNanull = zeros(geo.DofsPerElement,fe.GaussPointsPerElem,geo.NumElements);
                if this.fUseCrossFibreStiffness
                    a0a0n1 = anulldyadanull;
                    a0a0n2 = anulldyadanull;
                    a0base = anulldyadanull;
                end
                for m = 1 : geo.NumElements
                    u = geo.Nodes(:,geo.Elements(m,:));
                    for gp = 1 : fe.GaussPointsPerElem
                        
                        % Transform a0 to local fibre direction
                        pos = [0 fe.GaussPointsPerElem 2*fe.GaussPointsPerElem]+gp;
                        Jac = u*dNgp(:,pos);
                        if ismastercoord
                            loc_anull = Jac * anull(:,gp,m);
                        else
                            loc_anull = anull(:,gp,m);
                        end
                        loc_anull = loc_anull/norm(loc_anull);
                        
                        % Compute "the" two normals (any will do)
                        loc_anulln1 = circshift(loc_anull,1);
                        loc_anulln1 = loc_anulln1 - (loc_anulln1'*loc_anull) * loc_anull;
                        loc_anulln1 = loc_anulln1/norm(loc_anulln1);

                        loc_anulln2 = circshift(loc_anull,2);
                        loc_anulln2 = loc_anulln2 - (loc_anulln2'*loc_anull) * loc_anull - (loc_anulln2'*loc_anulln1) * loc_anulln1;
                        loc_anulln2 = loc_anulln2/norm(loc_anulln2);
                        
                        % forward transformation of a0 at gauss points
                        % (plotting only so far)
                        if ismastercoord
                            dNanull(:,gp,m) = dNgp(:,pos) * anull(:,gp,m);
                        else
                            dNanull(:,gp,m) = dNgp(:,pos) * (Jac \ anull(:,gp,m));
                        end
                        
                        % a0 dyad a0
                        pos = (m-1)*fe.GaussPointsPerElem+gp;
                        anulldyadanull(:,:,pos) = loc_anull*loc_anull';
                        a0base(:,:,pos) = [loc_anull loc_anulln1 loc_anulln2];
                        if this.fUseCrossFibreStiffness
                            a0a0n1(:,:,pos) = loc_anulln1*loc_anulln1';
                            a0a0n2(:,:,pos) = loc_anulln2*loc_anulln2';
                        end
                    end
                end
                this.a0oa0 = anulldyadanull;
                this.dNa0 = dNanull;
                this.a0Base = a0base;
                
                this.a0oa0n1 = [];
                this.a0oa0n2 = [];
                if this.fUseCrossFibreStiffness
                    this.a0oa0n1 = a0a0n1;
                    this.a0oa0n2 = a0a0n2;    
                end
            end
        end
        
        function initMuscleTendonRatios(this)
            mc = this.Model.Config;
            fe = mc.PosFE;
            this.HasTendons = ~isempty(mc.getTendonMuscleRatio(zeros(3,1)));
            tmr = zeros(fe.GaussPointsPerElem,fe.Geometry.NumElements);
            if this.HasTendons
                g = fe.Geometry;
                for m = 1:g.NumElements
                    % Get coordinates of gauss points in element
                    tmr(:,m) = mc.getTendonMuscleRatio(g.Nodes(:,g.Elements(m,:)) * fe.N(fe.GaussPoints));
                end
                this.MuscleTendonRatioNodes = mc.getTendonMuscleRatio(g.Nodes);
            end
            this.MuscleTendonRatioGP = tmr;
        end
    end
    
    methods(Access=private)
        function updateTendonMuscleParamsGP(this, mu)
            % Updates the markert coefficients for muscle or tendons
            % according to the current gauss points ratio of muscle and
            % tendon material.
            
            % All muscle - set without extra computations
            tmr = this.MuscleTendonRatioGP;
            if ~this.HasTendons
                tmr(:) = mu(9);
                this.MuscleTendonParamc10 = tmr;
                tmr(:) = mu(10);
                this.MuscleTendonParamc01 = tmr;
            else
                f = this.MuscleTendonParamMapFun;
                % Log-interpolated MR c10
                this.MuscleTendonParamc10 = f(tmr,mu(9),mu(11));
                % Log-interpolated MR c01
                this.MuscleTendonParamc01 = f(tmr,mu(10),mu(12));
            end
        end
    end
    
end
