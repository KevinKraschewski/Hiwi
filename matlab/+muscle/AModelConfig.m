classdef AModelConfig < handle
    %MODELCONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=private)
        Geometry;
        
        PosFE;
        
        PressFE;
    end
    
    methods
        function this = AModelConfig(geo)
            this.Geometry = geo;
            this.PosFE = triquadratic(geo);
            this.PressFE = trilinear(geo);
        end
    end
    
    methods(Access=protected)
        function anull = seta0(~, anull)
            % do nothing!
        end
        
        function x0 = getx0(~)
            x0 = [];
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(~, velo_dir, velo_dir_val)
            % set nothing!
        end
    end
    
    methods(Abstract, Access=protected)
        displ_dir = setPositionDirichletBC(this, displ_dir);
    end
    
    methods(Sealed)
        function [displ_dir, velo_dir, velo_dir_val] = getBC(this)
            N = this.PosFE.NumNodes;
            displ_dir = false(3,N);
            velo_dir = false(3,N);
            velo_dir_val = zeros(3,N);
            displ_dir = this.setPositionDirichletBC(displ_dir);
            [velo_dir, velo_dir_val] = this.setVelocityDirichletBC(velo_dir, velo_dir_val);
        end
        
        function anull = geta0(this)
            anull = zeros(3,this.Geometry.NumGaussp,this.PosFE.NumElems);
            anull = this.seta0(anull);
        end
    end
    
end
