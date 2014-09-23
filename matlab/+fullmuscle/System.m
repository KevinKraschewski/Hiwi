classdef System < muscle.System;
% System: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        num_motoneuron_dof;
        num_sarco_dof;
        num_spindle_dof;
        num_all_dof;
        off_moto;
        off_sarco;
        off_spindle;
        off_moto_full;
        off_sarco_full;
        off_spindle_full;
        
        input_motoneuron_link_idx;
        sarco_output_idx;
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
        
        noiseGen;
        
        % The upper limit polynomial for maximum mean current dependent on
        % the fibre type.
        %
        % Used at the assembly of B to provide suitable coefficient
        % functions.
        %
        % See also: models.motoneuron.experiments.ParamDomainDetection
        upperlimit_poly;
        
        % The offset for the sarcomere signals at t=0. Used to create
        % alpha(X,0) = 0
        %
        % See also: assembleX0
        sarco_mech_signal_offset
        
        Spindle;
        
        noise;
    end
    
    properties(Access=private)
        nfibres;
        basenoise;
    end
    
    methods
        function this = System(model)
            this = this@muscle.System(model);
            this.f = fullmuscle.Dynamics(this);
            
            % Load mean current limiting polynomial
            s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
            this.upperlimit_poly = s.upperlimit_poly;
            
            this.Inputs{1} = @(t)1;
            
            this.Spindle = fullmuscle.Spindle;
        end
        
        function configUpdated(this)
            ft = this.Model.Config.FibreTypes;
            nf = length(ft);
            this.nfibres = nf;
            
            configUpdated@muscle.System(this);
            
            %% Add input matrix B
            this.B = this.assembleB;
            
            %% Assemble noise signal for each fibre
            ng = models.motoneuron.NoiseGenerator;
            ng.setFibreType(ft(1));
            thenoise = zeros(nf,length(ng.indepNoise));
            thenoise(1,:) = ng.indepNoise;
            for k=2:nf
                ng.setFibreType(ft(k));
                thenoise(k,:) = ng.indepNoise;
            end
            this.noise = thenoise;
            this.basenoise = ng.baseNoise;
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.BaseDynSystem(this, mu, inputidx);
            
            if ~isempty(inputidx)
                % Create an input substitute that uses the true external
                % function and creates the effective noisy signal from it
                maxcurrents = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
                no = this.noise;
                bno = this.basenoise;
                
                % For nonconstant (=mu(4)) input use this
                %uin = this.Inputs{inputidx};
                
                ustr = '@(t)[bno(round(t)+1); ';
                for k=1:this.nfibres
                    % Precompute the factor as the input mean current will
                    % stay constant (at least from the external source)
                    uin = min(maxcurrents(k),mu(4));
                    rowfun = sprintf('no(%d,round(t)+1)*%g; ',k,uin);
                    
                    % For nonconstant (=mu(4)) input use this
                    %rowfun = sprintf('no(%d,round(t)+1)*min(maxcurrents(%d),uin(t))',k,k);
                    
                    ustr = [ustr rowfun];%#ok
                end
                ustr = [ustr ']'];
                this.u = eval(ustr);
            end
        end
        
        function uvwall = includeDirichletValues(this, t, uvw)
            uvwall_mech = includeDirichletValues@muscle.System(this, t, uvw(1:this.num_uvp_dof,:));
            uvwall = [uvwall_mech; uvw(this.num_uvp_dof+1:end,:)];
        end
        
        function [pm, h_mech] = plot(this, t, y, varargin)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('PM',[],@(v)isa(v,'PlotManager'));
            %i.addParamValue('F',[]);
            i.parse(varargin{:});
            r = i.Results;
            
            if ~isempty(r.PM)
                pm = r.PM;
            else
                pm = PlotManager(false,2,3);
                pm.LeaveOpen = true;
            end
            varargin(end+1:end+4) = {'Pool', false, 'PM', pm};
            [~, h_mech] = plot@muscle.System(this, t, y, varargin{:});
        end
        
        function v = coolExp(~, a, b, mu)
            v = exp(log(100)*mu)*(b-a)/100 + a;
        end
    end
    
    methods(Access=protected)
        
        function h = initRefinedPlot(this, t, y, r, pm)
            h = [];
            if length(t) > 1
                h = pm.nextPlot('signal','Motoneuron signal','t [ms]','V_m');
                pos = this.off_moto_full + (2:6:6*this.nfibres);
                vals = y(pos,:);
                axis(h,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h,'on');
                
                h2 = pm.nextPlot('force','Action potential','t [ms]','V_m');
                pos = this.off_sarco_full + (1:56:56*this.nfibres);
                vals = y(pos,:);
                axis(h2,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h2,'on');
                
                h3 = pm.nextPlot('force','Activation','t [ms]','A_s');
                pos = this.off_sarco_full + (53:56:56*this.nfibres);
                vals = y(pos,:);
                vals = bsxfun(@times, min(1,t), vals);
                axis(h3,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h3,'on');
                
                h4 = pm.nextPlot('spindle','Afferents','t [ms]','aff');
                pos = this.off_spindle_full + (1:9*this.nfibres);
                vals = this.Spindle.getAfferents(y(pos,:));
                axis(h4,[0 t(end) min(vals(:)) max(vals(:))]);
                hold(h4,'on');
                
                h5 = pm.nextPlot('spindle','Signal','t [ms]','aff');
                maxcurrents = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
                sp_sig = this.f.SpindleAffarentWeights * vals;
                sp_sig = min(maxcurrents - this.mu(4),sp_sig);
                axis(h5,[0 t(end) min(sp_sig(:)) max(sp_sig(:))]);
                hold(h5,'on');
                
                h = [h h2 h3 h4 h5];
            end
        end
        
        function refinedPlot(this, h, t, y, r, ts)
            if ts > 1
%                 cla(h(2:end));
                time_part = t(1:ts);
                pos = this.off_moto_full + (2:6:6*this.nfibres);
                signal = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
%                 cla(h(1));
                plot(h(1),time_part,signal);
                %plot(h,times,walpha,'LineWidth',2);
                
                pos = this.off_sarco_full + (1:56:56*this.nfibres);
                force = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
%                 cla(h(2));
                plot(h(2),time_part,force);
                
                pos = this.off_sarco_full + (53:56:56*this.nfibres);
                force = y(pos,1:ts);
                force = bsxfun(@plus, -this.sarco_mech_signal_offset, force);
                %force = bsxfun(@times, min(1,time_part), force);
                mc = this.Model.Config;
                walpha = mc.FibreTypeWeights(1,:,1) * force;
%                force = [force; walpha]
%                 cla(h(3));
                plot(h(3),time_part,force,'r',time_part,walpha,'b');
%                 plotyy(h(3),time_part,force,time_part,walpha);

                pos = this.off_spindle_full + (1:9*this.nfibres);
                affarents = this.Spindle.getAfferents(y(pos,1:ts));
                plot(h(4),time_part,affarents');
                
                maxcurrents = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
                sp_sig = this.f.SpindleAffarentWeights * affarents;
                sp_sig = min(maxcurrents - this.mu(4),sp_sig);
                plot(h(5),time_part,sp_sig);
            end
        end
        
        function updateDofNums(this, mc)
            updateDofNums@muscle.System(this, mc);
            
            this.num_motoneuron_dof = 6*this.nfibres;
            % Motoneurons are beginning after mechanics
            this.off_moto = this.num_uvp_dof; 
            this.off_moto_full = this.num_uvp_glob; 

            % Sarcomeres are beginning after motoneurons
            this.num_sarco_dof = 56*this.nfibres;
            this.off_sarco = this.off_moto + this.num_motoneuron_dof;
            this.off_sarco_full = this.off_moto_full + this.num_motoneuron_dof;
            
            % Spindles are beginning after sarcomeres
            this.num_spindle_dof = 9*this.nfibres;
            this.off_spindle = this.off_sarco + this.num_sarco_dof;
            this.off_spindle_full = this.off_sarco_full + this.num_sarco_dof;
            
            this.num_all_dof = this.off_spindle + this.num_spindle_dof;
            
            % Get the positions where the input signal is mapped to the
            % motoneurons
            this.input_motoneuron_link_idx = this.off_moto + (2:6:6*this.nfibres);
            
            this.sarco_output_idx = this.off_sarco + (53:56:56*this.nfibres);
        end
        
        function x0 = assembleX0(this)
            x0 = zeros(this.num_all_dof,1);
            % Get muscle x0
            x0(1:this.num_uvp_dof) = assembleX0@muscle.System(this);
            
            % Load dynamic/affine x0 coefficients for moto/sarco system
            % from file
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff.mat'));
            x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            
            ft = this.Model.Config.FibreTypes;
            smoff = zeros(this.nfibres,1);
            for k=1:this.nfibres
                x0ms = x0_motorunit.evaluate(ft(k));
                % add moto
                x0(this.off_moto + 6*(k-1) + (1:6)) = x0ms(1:6);
                % add sarco
                x0(this.off_sarco + 56*(k-1) + (1:56)) = x0ms(7:end);
                smoff(k) = x0ms(6+53);
                % add spindle
                x0(this.off_spindle + 9*(k-1) + (1:9)) = this.Spindle.y0;
            end
            this.sarco_mech_signal_offset = smoff;
        end
        
        function B = assembleB(this)
            % The divisor in both coefficient functions is the old para.CS
            % value!!
            i = this.input_motoneuron_link_idx;
            s = 1./(pi*this.coolExp(77.5e-4, 0.0113, this.Model.Config.FibreTypes).^2);

            B = sparse(i,ones(size(i)),s,this.num_all_dof,this.nfibres+1);
            for k=1:this.nfibres
                B(i(k),k+1) = s(k);%#ok
            end
            B = dscomponents.LinearInputConv(B);
%             B = dscomponents.AffLinInputConv;
%             % Base noise input mapping
%             B.addMatrix('1',sparse(i,ones(size(i)),s,this.num_all_dof,2));
%             
%             % Independent noise input mapping with µ_4 as mean current factor
%             % We need to restrict the maximum mean input current to a
%             % reasonable value for each fibre type. luckily, we know them
%             % already and hence can provide extra matrices with suitable
%             % coefficient functions.
%             %
%             % See also: 
%             maxvals = polyval(this.upperlimit_poly,this.Model.Config.FibreTypes);
%             for k=1:this.nfibres
%                 B.addMatrix(sprintf('min(%g,mu(4))',maxvals(k)),...
%                     sparse(i(k),2,s(k),this.num_all_dof,2));
%             end
        end
        
        function Daff = assembleDampingMatrix(this)
            Daff_mech = assembleDampingMatrix@muscle.System(this);
            
            Daff = dscomponents.AffLinCoreFun(this);
            extra = (6+56+9)*this.nfibres;
            D = blkdiag(Daff_mech.getMatrix(1),sparse(extra,extra));
            Daff.addMatrix('mu(1)',D);
        end
        
        function M = assembleMassMatrix(this)
            M = assembleMassMatrix@muscle.System(this);
            extra = (6+56+9)*this.nfibres;
            M = blkdiag(M,speye(extra));
        end
    end
    
end