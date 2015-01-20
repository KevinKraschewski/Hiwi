classdef MuscleTendonMix < muscle.AModelConfig
    % Muscle - Tendon mixed geometries
    %
    % Contains several examples:
    % Variant 1: Simple 
    
    properties
        Variant;
    end
    
    methods
        function this = MuscleTendonMix(variantnr)
            if nargin < 1
                variantnr = 1;
            end
            switch variantnr
                case 1
                    % Four cubes in a row
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,0:4,0:1);
                case {2 5}
                    % Single cube with same config as reference element
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,0:1,0:1);
                case {3 4}
                    % 2x4x2 geometry
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:2,0:2:6,0:2);
            end
            
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube27Node);
            this.Variant = variantnr;
            %this.NeumannCoordinateSystem = 'global';
        end
        
        function configureModel(this, m)
            m.T = 4;
            m.dt = .01;
            m.DefaultInput = 1;
            f = m.System.f;
            f.alpha = @(t)0;
            f.Pmax = 250;
            
            os = m.ODESolver;
            os.RelTol = .0001;
            os.AbsTol = .05;
            
            mu = [1; 0; 0; 0
                    %% Anisotropic parameters muscle+tendon (Markert)
                    4.02; 38.5; 7990; 16.6
                    %% Isotropic parameters muscle+tendon (Moonley-Rivlin)
                    35.6; 3.86; 2310; 1.15e-3]; % Micha
                    % 6.352e-10; 3.627 [Kpa] thomas alt
            switch this.Variant
                case {1 2}
                    mu(3) = 40;
                case 3
                    mu(3) = 10;
                case 4
                    mu(2) = m.T;
                    m.ODESolver.RelTol = .1;
                    m.ODESolver.AbsTol = .1;
                case 5
                    mu(2) = m.T;
                    m.ODESolver.RelTol = .1;
                    m.ODESolver.AbsTol = .1;
            end
            m.DefaultMu = mu;
            
            % Ramp up the external pressure
            m.System.Inputs{1} = this.getAlphaRamp(1,1);
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Overload this method to initialize model-specific quantities
            % that are fixed for each simulation
            %
            % Called by override of computeTrajectory in muscle.Model
            if any(this.Variant == [4 5])
                f = this.Model.System.f;
                f.alpha = this.getAlphaRamp(mu(2),1,0);
            end
        end
        
        function tmr = getTendonMuscleRatio(this)
            % Returns the [0,1] ratio between tendon and muscle at all
            % gauss points of all elements
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = getTendonMuscleRatio@muscle.AModelConfig(this);
            switch this.Variant
                case 1
                    tmr(:,1) = 0;%.4;
                    tmr(:,2) = .33;%.4;
                    tmr(:,3) = .66;%.4;
                    tmr(:,4) = 1;%.4;
                case {2 5} 
                    %tmr(:) = 1;
                    tmr(19:27,:) = 1;%.4;
                case {3 4}
                    f = @(x,z)min(1,max(0,-.5+.5./((x/2+.1).*(1.5*z+.08*x+.1))));
                    fges = @(x,z)f(x,z) + f(6-x,2-z);
                    fe = this.PosFE;
                    g = fe.Geometry;
                    for m = 1:g.NumElements
                        % Get coordinates of gauss points in element
                        gp = g.Nodes(:,g.Elements(m,:)) * fe.N(fe.GaussPoints);
                        tmr(:,m) = fges(gp(2,:),gp(3,:));
                    end
                    tmr = tmr/5;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            switch this.Variant
                % Pull on front face
                case {1 2 3}
                    if any(elemidx == [1 2 7 8]) && faceidx == 3
                        P = 1;
                    end
                %case 2
                %    if elemidx == 1 && faceidx == 1
                %        P = 1;
                %    end
            end
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            
            switch this.Variant
                case {1 2 5}
                    % Fix all on left and only the y,z directions of the back right
                    % nodes
                    displ_dir(2,geo.Elements(geo.NumElements,geo.MasterFaces(4,:))) = true;
                case 3
                    displ_dir(2,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                case 4
                    %displ_dir(2,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                    %displ_dir(2,geo.Elements([1 2 7 8],geo.MasterFaces(3,:))) = true;
                    displ_dir(:,geo.Elements([5 6 11 12],geo.MasterFaces(4,:))) = true;
                    displ_dir(:,geo.Elements([1 2 7 8],geo.MasterFaces(3,:))) = true;
                %case 5
                %    displ_dir(:,geo.Elements(geo.NumElements,geo.MasterFaces(4,:))) = true;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Variant
                case {1 2 5}
                    % Direction is y
                    anull(2,:,:) = 1;
                case {3 4}
                    %anull(2,:,:) = 1;
                    anull([2 3],:,:) = -1;
            end
        end
    end
    
end

