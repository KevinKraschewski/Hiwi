classdef triquadratic < fembase
    %TRILINEAR Summary of this class goes here
    %   Detailed explanation goes here
    %
    %% Cube node positions:
    % C1     -1    -1    -1
    % E2      0    -1    -1
    % C3      1    -1    -1
    % E4     -1     0    -1
    % E5      1     0    -1
    % C6     -1     1    -1
    % E7      0     1    -1
    % C8      1     1    -1
    % E9     -1    -1     0 
    % E10     1    -1     0
    % E11    -1     1     0
    % E12     1     1     0
    % C13    -1    -1     1
    % E14     0    -1     1
    % C15     1    -1     1
    % E16    -1     0     1
    % E17     1     0     1
    % C18    -1     1     1
    % E19     0     1     1
    % C20     1     1     1
    
    properties
        % cell array of dim n containing indices of all cubes that are
        % adjacent to the n-th point
        pts_cubes;
    end
    
    methods
        function this = triquadratic(geo)
            if nargin < 1
                geo = cubegeom;
            end
            this = this@fembase(geo);
            
            this.init;
        end
        
        function init(this)
            % Triquadratic basis functions
            %
            % N2 corner index
            % 1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20
            % 1   2     3   4             5     6        7     8  <-Corners
            % Combinatorial corner index
            % 1 2 3 4 6 7 8 9 10 12 16 18 19 20 21 22 24 25 26 27
            % (Missing 5,11,13,14,15,17,23)
            this.N = @(x)[(1-x(1,:)).*(1-x(2,:)).*(1-x(3,:)).*(-x(1,:)-x(2,:)-x(3,:)-2)/8;... % C1
                (1-x(1,:).^2).*(1-x(2,:)).*(1-x(3,:))/4;... % E2
                (1+x(1,:)).*(1-x(2,:)).*(1-x(3,:)).*(x(1,:)-x(2,:)-x(3,:)-2)/8;... % C3
                (1-x(2,:).^2).*(1-x(1,:)).*(1-x(3,:))/4;... % E4
                (1-x(2,:).^2).*(1+x(1,:)).*(1-x(3,:))/4;... % E5
                (1-x(1,:)).*(1+x(2,:)).*(1-x(3,:)).*(-x(1,:)+x(2,:)-x(3,:)-2)/8;... % C6
                (1-x(1,:).^2).*(1+x(2,:)).*(1-x(3,:))/4;... % E7
                (1+x(1,:)).*(1+x(2,:)).*(1-x(3,:)).*(x(1,:)+x(2,:)-x(3,:)-2)/8;... % C8
                (1-x(3,:).^2).*(1-x(1,:)).*(1-x(2,:))/4;... % E9
                (1-x(3,:).^2).*(1+x(1,:)).*(1-x(2,:))/4;... % E10
                
                (1-x(3,:).^2).*(1-x(1,:)).*(1+x(2,:))/4;... % E11
                (1-x(3,:).^2).*(1+x(1,:)).*(1+x(2,:))/4;... % E12
                (1-x(1,:)).*(1-x(2,:)).*(1+x(3,:)).*(-x(1,:)-x(2,:)+x(3,:)-2)/8;... % C13
                (1-x(1,:).^2).*(1-x(2,:)).*(1+x(3,:))/4;... % E14
                (1+x(1,:)).*(1-x(2,:)).*(1+x(3,:)).*(x(1,:)-x(2,:)+x(3,:)-2)/8;... %C15
                (1-x(2,:).^2).*(1-x(1,:)).*(1+x(3,:))/4;... % E16
                (1-x(2,:).^2).*(1+x(1,:)).*(1+x(3,:))/4;... % E17
                (1-x(1,:)).*(1+x(2,:)).*(1+x(3,:)).*(-x(1,:)+x(2,:)+x(3,:)-2)/8;... % C18
                (1-x(1,:).^2).*(1+x(2,:)).*(1+x(3,:))/4;... % E19
                (1+x(1,:)).*(1+x(2,:)).*(1+x(3,:)).*(x(1,:)+x(2,:)+x(3,:)-2)/8]; % C20
                
            this.gradN = @(x)[[-(1-x(2,:)).*(1-x(3,:)).*(-x(1,:)-x(2,:)-x(3,:)-2)-(1-x(1,:)).*(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(3,:)).*(-x(1,:)-x(2,:)-x(3,:)-2)-(1-x(1,:)).*(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(2,:)).*(-x(1,:)-x(2,:)-x(3,:)-2)-(1-x(1,:)).*(1-x(2,:)).*(1-x(3,:))]/8;... % 1
                [-2*x(1,:).*(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:).^2).*(1-x(3,:)) -(1-x(1,:).^2).*(1-x(2,:))]/4;...    
                [(1-x(2,:)).*(1-x(3,:)).*(x(1,:)-x(2,:)-x(3,:)-2)+(1+x(1,:)).*(1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(3,:)).*(x(1,:)-x(2,:)-x(3,:)-2)-(1+x(1,:)).*(1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(2,:)).*(x(1,:)-x(2,:)-x(3,:)-2)-(1+x(1,:)).*(1-x(2,:)).*(1-x(3,:))]/8;...
                [-(1-x(2,:).^2).*(1-x(3,:)) -2*x(2,:).*(1-x(1,:)).*(1-x(3,:))  -(1-x(2,:).^2).*(1-x(1,:))]/4;...
                [(1-x(2,:).^2).*(1-x(3,:)) -2*x(2,:).*(1+x(1,:)).*(1-x(3,:))  -(1-x(2,:).^2).*(1+x(1,:))]/4;...
                [-(1+x(2,:)).*(1-x(3,:)).*(-x(1,:)+x(2,:)-x(3,:)-2)-(1-x(1,:)).*(1+x(2,:)).*(1-x(3,:)) (1-x(1,:)).*(1-x(3,:)).*(-x(1,:)+x(2,:)-x(3,:)-2)+(1-x(1,:)).*(1+x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1+x(2,:)).*(-x(1,:)+x(2,:)-x(3,:)-2)-(1-x(1,:)).*(1+x(2,:)).*(1-x(3,:))]/8;...
                [-2*x(1,:).*(1+x(2,:)).*(1-x(3,:))  (1-x(1,:).^2).*(1-x(3,:)) -(1-x(1,:).^2).*(1+x(2,:))]/4;...
                [(1+x(2,:)).*(1-x(3,:)).*(x(1,:)+x(2,:)-x(3,:)-2)+(1+x(1,:)).*(1+x(2,:)).*(1-x(3,:)) (1+x(1,:)).*(1-x(3,:)).*(x(1,:)+x(2,:)-x(3,:)-2)+(1+x(1,:)).*(1+x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1+x(2,:)).*(x(1,:)+x(2,:)-x(3,:)-2)-(1+x(1,:)).*(1+x(2,:)).*(1-x(3,:))]/8;...
                [-(1-x(3,:).^2).*(1-x(2,:)) -(1-x(3,:).^2).*(1-x(1,:)) -2*x(3,:).*(1-x(1,:)).*(1-x(2,:))]/4;...
                [(1-x(3,:).^2).*(1-x(2,:)) -(1-x(3,:).^2).*(1+x(1,:)) -2*x(3,:).*(1+x(1,:)).*(1-x(2,:))]/4;... %E10
                
                [-(1-x(3,:).^2).*(1+x(2,:))  (1-x(3,:).^2).*(1-x(1,:)) -2*x(3,:).*(1-x(1,:)).*(1+x(2,:))]/4;...
                [(1-x(3,:).^2).*(1+x(2,:))  (1-x(3,:).^2).*(1+x(1,:)) -2*x(3,:).*(1+x(1,:)).*(1+x(2,:))]/4;...
                [-(1-x(2,:)).*(1+x(3,:)).*(-x(1,:)-x(2,:)+x(3,:)-2)-(1-x(1,:)).*(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:)).*(1+x(3,:)).*(-x(1,:)-x(2,:)+x(3,:)-2)-(1-x(1,:)).*(1-x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1-x(2,:)).*(-x(1,:)-x(2,:)+x(3,:)-2)+(1-x(1,:)).*(1-x(2,:)).*(1+x(3,:))]/8;...
                [-2*x(1,:).*(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:).^2).*(1+x(3,:))  (1-x(1,:).^2).*(1-x(2,:))]/4;...
                [(1-x(2,:)).*(1+x(3,:)).*(x(1,:)-x(2,:)+x(3,:)-2)+(1+x(1,:)).*(1-x(2,:)).*(1+x(3,:)) -(1+x(1,:)).*(1+x(3,:)).*(x(1,:)-x(2,:)+x(3,:)-2)-(1+x(1,:)).*(1-x(2,:)).*(1+x(3,:)) (1+x(1,:)).*(1-x(2,:)).*(x(1,:)-x(2,:)+x(3,:)-2)+(1+x(1,:)).*(1-x(2,:)).*(1+x(3,:))]/8;...
                [-(1-x(2,:).^2).*(1+x(3,:)) -2*x(2,:).*(1-x(1,:)).*(1+x(3,:))   (1-x(2,:).^2).*(1-x(1,:))]/4;...
                [(1-x(2,:).^2).*(1+x(3,:)) -2*x(2,:).*(1+x(1,:)).*(1+x(3,:))   (1-x(2,:).^2).*(1+x(1,:))]/4;...
                [-(1+x(2,:)).*(1+x(3,:)).*(-x(1,:)+x(2,:)+x(3,:)-2)-(1-x(1,:)).*(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(3,:)).*(-x(1,:)+x(2,:)+x(3,:)-2)+(1-x(1,:)).*(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(2,:)).*(-x(1,:)+x(2,:)+x(3,:)-2)+(1-x(1,:)).*(1+x(2,:)).*(1+x(3,:))]/8;...
                [-2*x(1,:).*(1+x(2,:)).*(1+x(3,:))  (1-x(1,:).^2).*(1+x(3,:))  (1-x(1,:).^2).*(1+x(2,:))]/4;...
                [(1+x(2,:)).*(1+x(3,:)).*(x(1,:)+x(2,:)+x(3,:)-2)+(1+x(1,:)).*(1+x(2,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(3,:)).*(x(1,:)+x(2,:)+x(3,:)-2)+(1+x(1,:)).*(1+x(2,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(2,:)).*(x(1,:)+x(2,:)+x(3,:)-2)+(1+x(1,:)).*(1+x(2,:)).*(1+x(3,:))]/8];
            
            % Transform cube data to dof/element data
            p = this.geo.pts;
            c = this.geo.cubes;
            nc = size(c,1);
            % Transformation matrix for corners to corner+edges locations
            i = [1 2  2  3  4  4  5  5 6  7  7 8  9  9 10 10 11 11 12 12 13 14 14 15 16 16 17 17 18 19 19 20];  
            j = [1 1  2  2  1  3  2  4 3  3  4 4  1  5 2  6  3  7  4  8  5  5  6  6  5  7  6  8  7  7  8  8];
            s = [1 .5 .5 1 .5 .5 .5 .5 1 .5 .5 1 .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 1  .5 .5 .5 .5 1  .5 .5 1];
            T = sparse(i,j,s,20,8);
            nodes = zeros(3,nc*20);
            elems = zeros(nc,20);
            % Iterate all cubes and collect nodes
            for cidx = 1:nc
                cube = c(cidx,:);
                pos = 20*(cidx-1)+1:20*cidx;
                nodes(:,pos) = p(:,cube)*T';
                elems(cidx,:) = pos;
            end
            [nodes, ~, elems] = unique(nodes','rows');
            this.nodes = nodes';
            this.elems = reshape(elems,20,[])';
            
            init@fembase(this);
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(this.elems,1)
                hlp = this.elems(i,[1 2 1 4 1 9 2 3 3 5 3 10 4 6 5 8 6 7 ...
                    6 11 7 8 8 12 9 13 10 15 11 18 12 20 13 14 13 16 ...
                    14 15 15 17 16 18 17 20 18 19 19 20]);
                e(end+1:end+24,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows');
            this.edges = e;
        end
    end
    
    methods(Static)
        function res = test_QuadraticBasisFun
            q = triquadratic;
            res = fembase.test_BasisFun(q);
            
            % test for correct basis function values on nodes
            [X,Y,Z] = ndgrid(-1:1:1,-1:1:1,-1:1:1);
            p = [X(:) Y(:) Z(:)]';
            % Remove 7 (inner points not used here)
            p(:,[5,11,13,14,15,17,23]) = [];
            res = res && isequal(q.N(p),eye(20));
        end
    end
    
end

