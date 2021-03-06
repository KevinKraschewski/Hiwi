classdef Dynamics < muscle.Dynamics;
% Dynamics: 
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
    
    %% Properties for motoneuron - sarcomere linking
    properties
        % The `V_s` value of the motoneuron input at which the MSLink_MaxFactor should be attained
        %
        % @type double
        MSLink_MaxFactorSignal = 40;
        
        % The maximal factor with which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations
        %
        % @type double
        MSLink_MaxFactor = 7;
        
        % The minimal factor at which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations.
        %
        % Between both limits, an exponential Gaussian weight is applied with a certain radius,
        % so that the amplification factor has its minimum value at zero and maximum value at
        % MSLink_MaxFactorSignal.
        %
        % See also FibreDynamics.getLinkFactor
        %
        % @type double
        MSLink_MinFactor = .3;
        
        MSLinkFun;
        MSLinkFunDeriv;
        
        moto_sarco_link_moto_out;
        moto_sarco_link_sarco_in;
        spindle_moto_link_moto_in;
        
        FrequencyDetector;
        
        % The factors with which the primary and secondary affarent of the
        % spindle is multiplied before considered a "mean input current"
        % (then added to the external signal)
        %
        % @type rowvec<double> @default [1 1]*0.002
        SpindleAffarentWeights = sparse([1 1]*0.002);
    end
    
    properties(Dependent)
        UseFrequencyDetector;
    end

    properties(Access=private)
        % The upper limit for the mean input current fed to the motoneuron
        % soma. as in this current version this is the sum of spindle
        % feedback and external signal, the max value needs to be available
        % here to limit the sum.
        % The external signal is added at an higher level in the ode (B*u
        % component), but yet the sum of (spindle+ext_sig) <
        % max_moto_signal, which is why the external signal is accessed
        % here, too.
        max_moto_signals;
        
        fUseFD = false;
        freq_kexp;
        sarcoconst1;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@muscle.Dynamics(sys);
            s = load(fullfile(fileparts(which('fullmuscle.Model')),'FrequencyKexp.mat'));
            this.freq_kexp = s.kexp.toTranslateBase;
        end
        
        function configUpdated(this)
            sys = this.System;
            mc = sys.Model.Config;
            ft = mc.FibreTypes;
            this.nfibres = length(ft);
            
            this.lambda_dot_pos = double.empty(2,0);
            this.FrequencyDetector = fullmuscle.FrequencyDetector(this.nfibres);
            if sys.HasSpindle
                this.lambda_dot_pos = mc.SpindlePositions;
                this.lambda_dot = zeros(1,this.nfibres);
            end
            
            % Needed in evaluate & getStateJacobian
            this.sarcoconst1 = sys.Sarcomere.SarcoConst(1,:);
            
            configUpdated@muscle.Dynamics(this);
            
            % fDim and xDim are from muscle.Dynamics, so add moto+sarco
            this.fDim = this.fDim + (6+56)*this.nfibres;
            this.xDim = this.xDim + (6+56)*this.nfibres;
            if sys.HasSpindle
                this.fDim = this.fDim + 9*this.nfibres;
                this.xDim = this.xDim + 9*this.nfibres;
            end
            
            this.moto_sarco_link_moto_out = sys.off_moto + (2:6:6*this.nfibres);
            this.moto_sarco_link_sarco_in = sys.off_sarco + (1:56:56*this.nfibres);
            this.spindle_moto_link_moto_in = this.moto_sarco_link_moto_out;
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@muscle.Dynamics(this, mu);
            % Create function handles for the link function
            diff = this.MSLink_MaxFactor-this.MSLink_MinFactor;
            funstr = sprintf('%g + exp(-(x-%g).^2/150)*%g',...
                this.MSLink_MinFactor,...
                this.MSLink_MaxFactorSignal,diff);
            this.MSLinkFun = eval(['@(x)' funstr ';']);
            funstr = sprintf('-exp(-(x-%g).^2/150).*(%g*(x-%g))/75',...
                this.MSLink_MaxFactorSignal,...
                diff,...
                this.MSLink_MaxFactorSignal);
            this.MSLinkFunDeriv = eval(['@(x)' funstr ';']);
            
            sys = this.System;
            % Register the ODE callback for the frequency integration
            slv = sys.Model.ODESolver;
            if ~isa(slv,'solvers.MLWrapper')
                error('Only programmed to work with ML-builtin solvers so far!');
            end
            if this.fUseFD
                fd = this.FrequencyDetector;
                slv.odeopts = odeset(slv.odeopts,...
                    'OutputFcn',@(t,y,flag)fd.processSignal(t,y'),...
                    'OutputSel',this.moto_sarco_link_moto_out);
            else
                slv.odeopts.OutputFcn = [];
            end
            this.max_moto_signals = polyval(sys.upperlimit_poly,sys.Model.Config.FibreTypes);
        end
        
        function dy = evaluate(this, y, t)
            sys = this.System;
            dy = zeros(this.fDim,1);
            
            %% Mechanics
            uvp_pos = 1:sys.NumTotalDofs;
            % Use uvp as argument and also pass in s (=sarco forces)
            force = max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset);
%             force = y(sys.sarco_output_idx);
%             force = zeros(length(sys.sarco_output_idx),1);
            uvps = [y(uvp_pos); force];
            dy(uvp_pos) = evaluate@muscle.Dynamics(this, uvps, t);
            
            %% Motoneurons
            mo = sys.Motoneuron;
            moto_pos = sys.off_moto+(1:sys.num_motoneuron_dof);
            dy_m = mo.dydt(reshape(y(moto_pos),6,[]),t);
            dy(moto_pos) = dy_m(:);
            
            %% Sacromeres
            sa = sys.Sarcomere;
            sarco_pos = sys.off_sarco + (1:sys.num_sarco_dof);
            ys = reshape(y(sarco_pos),56,[]);
            dys = sa.dydt(ys, t);
            dy(sarco_pos) = dys(:);
            
            %% Link of motoneurons to sarcomeres
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            signal = fac.*moto_out./this.sarcoconst1';
            % Add signal to corresponding locations
            dy(this.moto_sarco_link_sarco_in) = dy(this.moto_sarco_link_sarco_in) + signal;
            
            if sys.HasSpindle
                %% Spindles
                sp = sys.Spindle;
                spindle_pos = sys.off_spindle + (1:sys.num_spindle_dof);
                yspindle = reshape(y(spindle_pos),9,[]);

                % Link of spindle to motoneuron
                % Get single spindle signals
                spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                % Compute the mean value over all signals
                spindle_sig = ones(1,this.nfibres)*mean(spindle_sig);
                % Get current external signal
                ext_sig = sys.Inputs{2,1}(t);
                % Use the upper bounded sum
                eff_spindle_sig = min(spindle_sig,this.max_moto_signals - ext_sig);
                % Compute noisy signal
                noise_sig = mo.TypeNoise(:,round(t)+1)'.*eff_spindle_sig.*mo.FibreTypeNoiseFactors;
    %             fprintf('Spindle->Neuron: adding %g at dy(%d)\n',noise_sig,this.spindle_moto_link_moto_in);
                dy(this.spindle_moto_link_moto_in) = ...
                    dy(this.spindle_moto_link_moto_in) + noise_sig';

                % Spindle actual
                % Get motoneuron frequency
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else
                    freq = this.freq_kexp.evaluate([sys.Model.Config.FibreTypes; eff_spindle_sig+ext_sig]);
                end
                dys = sp.dydt(yspindle,t,freq,this.lambda_dot,0);
                dy(spindle_pos) = dys(:);
            end
        end
        
        function [J, JLamDot] = getStateJacobian(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;
            sys = this.System;
            
            %% Mechanics
            uvp_pos = 1:sys.NumTotalDofs;
            force = max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset);
%             force = y(sys.sarco_output_idx);
%             force = zeros(length(sys.sarco_output_idx),1);
            uvps = [y(uvp_pos); force];
            [J, Jalpha, JLamDot] = getStateJacobian@muscle.Dynamics(this, uvps, t);
            
            %% Motoneuron
            mo = sys.Motoneuron;
            for k=1:this.nfibres
                moto_pos = sys.off_moto + 6*(k-1) + (1:6);
                J = blkdiag(J,mo.Jdydt(y(moto_pos),t,k));
            end
            
            %% Sarcomeres
            sa = sys.Sarcomere;
            for k=1:this.nfibres
                sarco_pos = sys.off_sarco + 56*(k-1) + (1:56);
                J = blkdiag(J,sa.Jdydt(y(sarco_pos),t,k));
            end
            
            %% Motoneuron to Sarcomere coupling
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            dfac = this.MSLinkFunDeriv(moto_out);
            dsignal_dmotoout = (dfac .* moto_out + fac)./this.sarcoconst1';
            for k=1:this.nfibres
                J(this.moto_sarco_link_sarco_in(k),this.moto_sarco_link_moto_out(k)) = dsignal_dmotoout(k);
            end
            
            %% Sarcomere to mechanics coupling
            % The JS matrix is generated during the computation of the
            % mechanics jacobian, as the element/gauss loop is computed
            % there anyways. its not 100% clean OOP, but works for now.
            J(sys.NumStateDofs + (1:sys.NumDerivativeDofs), sys.off_sarco+(1:sys.num_sarco_dof)) = Jalpha;
            
            %% Spindle stuff
            if sys.HasSpindle
                sp = sys.Spindle;
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else
                    % For no detection, the current spindle signal is
                    % required
                    spindle_pos = sys.off_spindle + (1:sys.num_spindle_dof);
                    yspindle = reshape(y(spindle_pos),9,[]);
                    % Get single spindle signals
                    spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                    % Compute the mean value over all signals
                    spindle_sig = ones(1,this.nfibres)*mean(spindle_sig);
                    % Use the upper bounded sum
                    eff_spindle_sig = min(this.max_moto_signals, spindle_sig + sys.Inputs{2,1}(t));
                    freq_kexp_arg = [sys.Model.Config.FibreTypes; eff_spindle_sig];
                    freq = this.freq_kexp.evaluate(freq_kexp_arg);
                end
                
                i = []; j = []; s = [];
                moto_pos = 2:6:6*this.nfibres;
                for k=1:this.nfibres
                    spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                    
                    %% Spindles by themselves
                    [Jspin, Jspin_dLdot, Jspin_dmoto] = sp.Jdydt(y(spindle_pos), t, freq(k), this.lambda_dot(k), 0);
                    J = blkdiag(J,Jspin);
                    
                    %% Mechanics to spindle coupling
                    J(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = Jspin_dLdot'*JLamDot(k,:);
                    
                    %% Spindle to Motoneuron coupling
                    daffk_dy = this.SpindleAffarentWeights*sp.getAfferentsJacobian(y(spindle_pos));
                    dnoise_daff = mo.TypeNoise(:,round(t)+1).*mo.FibreTypeNoiseFactors(:);
                    i = [i repmat(moto_pos',1,9)];%#ok
                    j = [j repmat(9*(k-1) + (1:9),this.nfibres,1)];%#ok
                    s = [s dnoise_daff*daffk_dy/this.nfibres];%#ok
                    
                    %% Moto to Spindle coupling for learned frequencies
                    if ~this.fUseFD
                        kexp_Jac = this.freq_kexp.getStateJacobian(freq_kexp_arg(:,k));
                        
                        J(spindle_pos,spindle_pos) = J(spindle_pos,spindle_pos) ...
                            + Jspin_dmoto'*kexp_Jac(2)*daffk_dy;
                    end
                end
                J(sys.off_moto + (1:sys.num_motoneuron_dof),...
                  sys.off_spindle + (1:sys.num_spindle_dof))...
                    = sparse(i(:),j(:),s(:),sys.num_motoneuron_dof,sys.num_spindle_dof);
            end
        end
        
        function res = test_Jacobian(this, y, t, mu)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            %
            % Currently the tests using viscosity are commented out as we
            % assume linear damping, which is extracted as extra `A(t,\mu)`
            % part in the models' system
            sys = this.System;
            if nargin < 4
                mu = sys.Model.DefaultMu;
                if nargin < 3
                    t = 1000;
                    if nargin < 2
                        y = this.System.x0.evaluate(mu);
                    end
                end
            end
            
            % Also test correct computation of JLamDot
            d = size(y,1);
            dx = ones(d,1)*sqrt(eps(class(y))).*max(abs(y),1);
            sys.prepareSimulation(mu, sys.Model.DefaultInput);
            this.evaluate(y,t);
            ldotbase = this.lambda_dot;
            uv = sys.NumStateDofs+sys.NumDerivativeDofs;
            LAM = repmat(ldotbase',1,uv);
            dlam = zeros(this.nfibres,uv);
            for k = 1:uv
                h = zeros(d,1);
                h(k) = dx(k);
                this.evaluate(y+h,t);
                dlam(:,k) = this.lambda_dot;
            end
            JLamFD = (dlam - LAM)./dx(1:uv)';
            [~, JLamDot] = this.getStateJacobian(y,t);
            difn = norm(JLamFD-JLamDot);
            res = difn < 1e-13;
            
            freq = ones(1,this.nfibres)*30;
            dx = ones(this.nfibres,1)*sqrt(eps(class(ldotbase))).*max(abs(ldotbase),1);
            sp = sys.Spindle;
            for k = 1:this.nfibres
                spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                fx = sp.dydt(y(spindle_pos),t,freq(k),ldotbase(k),0);
                fxh_dldot = sp.dydt(y(spindle_pos),t,freq(k),ldotbase(k)+dx(k),0);
                fxh_dmoto = sp.dydt(y(spindle_pos),t,freq(k)+dx(k),ldotbase(k),0);
                Jspin_Ldot_FD = (fxh_dldot-fx) / dx(k);
                Jspin_moto_FD = (fxh_dmoto-fx) / dx(k);
                [~, Jspin_dLdot, Jspin_dmoto] = sp.Jdydt(y(spindle_pos), t, freq(k), ldotbase(k), 0);
                difn = norm(Jspin_Ldot_FD - Jspin_dLdot');
                res = res && difn < 1e-7;
                difn = norm(Jspin_moto_FD - Jspin_dmoto');
                res = res && difn < 1e-7;
            end
            
            res = res & test_Jacobian@muscle.Dynamics(this, y, t, mu);
        end
    end
    
    %% Getter / Setter
    methods
        function value = get.UseFrequencyDetector(this)
            value = this.fUseFD;
        end
        
        function set.UseFrequencyDetector(this, value)
            if ~isequal(value,this.fUseFD) || isempty(this.fUseFD)
                this.fUseFD = value;
                % Update pattern (learned frequency causes different
                % spindle pattern)
                if ~isempty(this.System.Model.Config)
                    this.JSparsityPattern = this.computeSparsityPattern;
                end
            end
        end
    end
    
    methods(Access=protected)
        function SP = computeSparsityPattern(this)
            [SP, SPalpha, SPLamDot] = computeSparsityPattern@muscle.Dynamics(this);
            sys = this.System;
            
            % Neuro
            i = [1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6];
            j = [1,2,1,2,3,4,5,6,2,3,2,4,2,5,2,6];
            J_moto = sparse(i,j,true,6,6);
            for k=1:this.nfibres
                SP = blkdiag(SP, J_moto);
            end
            
            % Sarco
            J_sarco = sys.Sarcomere.JSparsityPattern;
            for k=1:this.nfibres
                SP = blkdiag(SP, J_sarco);
            end
            
            if sys.HasSpindle
                % Spindle
                JSp = sys.Spindle.JSparsityPattern;
                for k=1:this.nfibres
                    SP = blkdiag(SP, JSp);
                end
            end
            
            % Moto -> Sarco link
            for k=1:this.nfibres
                % first entry of sarco gets 2nd output of motoneuron
                off_sarco = sys.off_sarco + (k-1)*56 + 1;
                off_moto = sys.off_moto + (k-1)*6 + 2;
                SP(off_sarco,off_moto) = true;
            end
                    
            % Sarco -> Mechanics
            SP(sys.NumStateDofs + (1:sys.NumDerivativeDofs), sys.off_sarco+(1:sys.num_sarco_dof)) = SPalpha;
            
            if sys.HasSpindle
                % Spindle -> Motoneuron link
                i = []; j = [];
                moto_pos = 2:6:6*this.nfibres;
                for k=1:this.nfibres
                    i = [i repmat(moto_pos',1,9)];%#ok
                    j = [j repmat(9*(k-1) + (1:9),this.nfibres,1)];%#ok
                end
                SP(sys.off_moto + (1:sys.num_motoneuron_dof),...
                  sys.off_spindle + (1:sys.num_spindle_dof))...
                    = sparse(i(:),j(:),ones(numel(i),1),...
                    sys.num_motoneuron_dof,sys.num_spindle_dof);
                
                % Mechanics -> spindle
                Jspin_Ldot = double(sys.Spindle.JLdotSparsityPattern);
                Jspin_dmoto = double(sys.Spindle.JMotoSparsityPattern);
                Jspin_Aff = double(sys.Spindle.JAfferentSparsityPattern);
                for k=1:this.nfibres
                    spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                    SP(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = ...
                        logical(Jspin_Ldot*double(SPLamDot(k,:)));
                    % Create connecting link to self only when no frequency
                    % detector is used!
                    if ~this.fUseFD
                        SP(spindle_pos,spindle_pos) = ...
                            SP(spindle_pos,spindle_pos) | logical(Jspin_dmoto*any(Jspin_Aff));
                    end
                end
            end
        end
    end
    
    methods(Static)
        function res = test_Dynamics
            m = fullmuscle.Model(fullmuscle.CPull(1));
            f = m.System.f;
            res = true;
            
            f.UseFrequencyDetector = false;
            [t, y] = m.simulate;
            res = res & f.test_Jacobian;
            res = res & f.test_Jacobian(y(:,end),t(end));
            
            f.UseFrequencyDetector = true;
            [t, yfd] = m.simulate;
            res = res & f.test_Jacobian;
            res = res & f.test_Jacobian(yfd(:,end),t(end));
        end
    end
end