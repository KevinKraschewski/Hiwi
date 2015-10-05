classdef ShakerDefault < models.muscle.AMuscleConfig
    
    
    properties(SetAccess=private)
        ylen;
        radfun;
        stretchfun;
        %New Parameters for frequency, a tolerance for the maxDy of
        %TMRFunction (for the discretisation part) and the maximal amount
        %of points that are usable
        frequency;
        TOL;
        maxYPoints; %Default als 0! Fuer die Disc Methode (Abfrage ob der Wert 0 ist)
    end
    
    methods
        function this = ShakerDefault(varargin)
            this = this@models.muscle.AMuscleConfig(varargin{:});
           
            
            %Default Stretch is exp
            this.addOption('Stretch','exp');
            %Default frequency is 50Hz
            this.addOption('Frequency',50);
            %Default tolerance 0.1
            this.addOption('TOL',10^(-1));
            %Default maxYPoints 0, as this value is used in a method!
            this.addOption('maxYPoints',0);
            
            
            
            % Compute outer shape
            this.ylen = 100;
            k = kernels.GaussKernel(20);
            this.radfun = @(x)10*k.evaluate(this.ylen/2,x)+2;
            
            % Also invokes geo setup
            this.init;
            
            %Sets the frequency property
            this.frequency=this.Options.Frequency;
            %Sets the tolerance property
            this.TOL=this.Options.TOL;
            %Sets the maxYPoints property
            this.maxYPoints=this.Options.maxYPoints;
            
            ref = [1.03, 1.4];
            %Man kann nun ueber 'Gauss Gamma' einen GaußKernel als
            %Stretchfunktion verwenden, wichtig ist das Leerzeichen
            %zwischen Gauss und dem gewuenschten Gamma.
            %Bei groesseren Gamma -> Ungenauigkeiten im maximal Wert, fuehrt zur Unbrauchbarkeit!
            if length(this.Options.Stretch) >= 5 && strcmp(this.Options.Stretch(1:5),'Gauss')
                Gamma=str2double(this.Options.Stretch(6:length(this.Options.Stretch)));
                k=kernels.GaussKernel(Gamma);
                kexp=kernels.KernelExpansion;
                kexp.Kernel=k;
                kexp.Centers.xi=1;       % Funktion soll auf [0,1] definiert sein
                kexp.Ma=(ref(2)-ref(1)); %Damit insgesamt das Maximum wieder bei ref(2) liegt
               stretchfun_0_1=@(x) (ref(1)-kexp.evaluate(0))+kexp.evaluate(x); %Der erste Term ist uA da um "Ungenauigkeiten" auszugleichen
            else
            switch this.Options.Stretch
                case 'exp'
                    stretchfun_0_1 =@(x) ref(1)+(ref(2)-ref(1)).*x.^4;
                case 'log'
                    stretchfun_0_1 =@(x)((ref(2)-ref(1))/log(51)).*log(50.*(1/50+x))+ref(1);
                case 'linear'
                    stretchfun_0_1 =@(x) (ref(2)-ref(1)).*x+ref(1);
            end
            end
            len = this.ylen;
            % Transfer the [0,1] argument version into the given domain
            % length (mirror at half)
            this.stretchfun = @(y)(y<=len/2) .* (stretchfun_0_1(2*y/len)) ...
                + (y>len/2) .* (stretchfun_0_1(1-2*(y-len/2)/len));
            
            this.VelocityBCTimeFun = general.functions.Sinus(this.frequency); %gewuenschte Frequenz
        end
        
        function plotTMR(this)
            y = linspace(0,this.ylen,1000);
            tmr = this.getTendonMuscleRatio([zeros(size(y)); y]);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('tmr','Prescribed stretches and resulting TMR','y-axis [mm]','Value');
            plot(ax, y,this.stretchfun(y),'r',y,tmr,'b');
            legend('Stretch','TMR');
            pm.done;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 75;
            m.dt = .5;
            
            mu = m.DefaultMu;
            % Small viscosity
            mu(1) = .005;
            m.DefaultMu = mu;
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            
            f = this.Model.System.f;
            % This is empty if no prepareSimulation has been called yet
            if isempty(f.AnisoPassiveTendon)
                tmr = zeros(1,size(points,2));
                return;
            end
            
            % Berechne zunaechst die Kraft die an der Sehne aufgebracht
            % werden muss, um die bei x=0 vorgegebene Sehnenverzerrung zu erreichen
            
            F = f.AnisoPassiveTendon(this.stretchfun(0))*pi*this.radfun(0)^2;
            
            % Fläche als Funktion des Ortes
            A = @(y)this.radfun(y).^2*pi;
            
            % Extrahiere y-Koordinaten aus den gegebenen Punkten
            ycoord = points(2,:);
            
            % Berechne Stress der einzelkomponenten bei gegebenem Stretch
            % an den gegebenen x-Positionen
            musclestress = f.AnisoPassiveMuscle(this.stretchfun(ycoord));
            tendonstress = f.AnisoPassiveTendon(this.stretchfun(ycoord));
            
            % Umstellen von
            % F = \sigma*A = (1-r)\sigma_m + r\sigma_t
            % nach r ergibt
            % r = (F/A - }\sigma_m) / (sigma_t - sigma_m)
            tmr = (F./A(ycoord) - musclestress)./(tendonstress-musclestress);
        end
        
        function [ypoints] = getDiscr(this,TMRFunc) %TMRFunc muss hier eine Funktion sein!!
            this.ylen=this.ylen/2;
            this.maxYPoints=160; %Nur zum Testen!
            this.maxYPoints=this.maxYPoints/2
            this.TOL=10^(-2); %Nur zum Testen
            deltaX=this.ylen/4; %Anfangsschrittweite           
            maxdy=this.TOL; %Default maximales dy als Differenzen der TMRFunktionswerte

            
            if this.maxYPoints~=0 %Anzahl an Intervallpunkten vorgegeben
                ypoints=zeros(this.maxYPoints,1); %Initialisieren
                %%Schleife ueber die vorgegebene Anzahl an Punkten
                %%abzueglich dem letzten Punkt, da dieser Fest ist
                for i=1 : this.maxYPoints-1 %Erster und letzter Punkt fest
                    
                    %Solange man mit dem naechsten Punkt rechts rausfallen
                    %wuerde
                    while ypoints(i)+deltaX > this.ylen
                        deltaX=deltaX/2;
                    end
                    
                    %%Schleife bis Tolerierte Schrittweite erreicht wird
                    while 1 %Solange die tolerierte Schrittweite nicht erreicht wird repeat
                        if abs(TMRFunc(ypoints(i))-TMRFunc(ypoints(i)+deltaX))<this.TOL;%Schrittweite akzeptiert
                            ypoints(i+1)=ypoints(i)+deltaX;                        %Naechster Punkt ist der alte plus ermitteltes "notwendiges" deltaT 
                            break;
                        else %Schrittweite nicht akzeptiert -> Versuche es mit halber Schrittweite
                       deltaX=deltaX/2;
                        end                        
                    end
                    deltaX=this.ylen/4; %Schrittweite am Ende immer auf Anfangsschrittweite zuruecksetzen
                    %Letzte Abfrage ob insgesamt die Toleranz erreicht
                    %wurde
                    if abs(TMRFunc(this.ylen)-TMRFunc(ypoints(i)))<this.TOL; %Irgendwann zuvor Toleranz erreicht -> fertig
                       maxdy=abs(TMRFunc(this.ylen)-TMRFunc(ypoints(i)));
                       ypoints(i+1)=this.ylen;
                       ypoints=ypoints(1:i+1); %Damit nicht x mal der selbe Punkt verwendet wird
                       disp(['The Tolerance was already reached with the number of ' num2str(i+1) ' points and the max of dy is ' num2str(maxdy)]);
                       help=0;
                       break
                    end                    
                end
                if help && i==this.maxYPoints-1 && abs(TMRFunc(this.ylen)-TMRFunc(ypoints(this.maxYPoints-1)))< this.TOL; %Punkte reichen genau aus
                maxdy=abs(TMRFunc(this.ylen)-TMRFunc(ypoints(this.maxYPoints-1)));
                ypoints(this.maxYPoints)=this.ylen;
                disp(['The Tolerance was reached with the amount of ' num2str(this.maxYPoints) ' points and the max of dy is ' num2str(maxdy)]);
                
                elseif i==this.maxYPoints-1 && abs(TMRFunc(this.ylen)-TMRFunc(ypoints(this.maxYPoints-1))) > this.TOL; %Nicht genug Punkte fuer gewuenschte Toleranz
                maxdy=abs(TMRFunc(this.ylen)-TMRFunc(ypoints(this.maxYPoints-1)));
                disp(['The Tolerance can not be reached with this amount of points. The max of dy would be ' num2str(maxdy)]); %Error oder disp?
                ypoints(i+1)=this.ylen;
                end
                
            else %Nur Toleranz vorgegeben -> So viele Punkte wie "noetig".
                ypoints=zeros(1000,1); %Preallocation for better performance
                i=1;
                this.TOL = 10^(-2); %Nur zum Testen!
                while 1 %Endless loop until tolerance globally reached
                    
                    %Solange man mit dem naechsten Punkt rechts rausfallen
                    %wuerde
                    while ypoints(i)+deltaX > this.ylen
                        deltaX=deltaX/2;
                    end                    
                    if (abs(TMRFunc(ypoints(i))-TMRFunc(this.ylen))) < this.TOL %Toleranz insgesamt erreicht -> Break
                        ypoints(i+1)=this.ylen; %Letzter Punkt fix
                        ypoints=ypoints(1:i+1); %Nur die wirklich gebrauchten Punkte verwenden, restl 0en der Preallocation "wegschneiden"
                        disp(['Tolerance reached and used ' num2str(i) 'points, maxDy is ' num2str(maxdy)]);
                        break
                    else
                        while abs(TMRFunc(ypoints(i))-TMRFunc(ypoints(i)+deltaX)) > this.TOL %Verkleinere Schrittweite bis passt
                        deltaX=deltaX/2;
                        end
                        maxdy=abs(TMRFunc(ypoints(i))-TMRFunc(ypoints(i)+deltaX));
                        ypoints(i+1)=ypoints(i)+deltaX;
                    end
                    i=i+1;
                end
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
%             m  = models.muscle.Model(models.muscle.examples.Belly);
%             mu = m.DefaultMu;
%             
%             S_m = models.muscle.functions.MarkertLawOriginal(mu(5),mu(6));
%             AnisoPassiveMuscle = S_m.getFunction;
%             S_t = general.functions.CubicToLinear(mu(7),mu(8));
%             AnisoPassiveTendon = S_t.getFunction;
%             %Stretchfun und radfun fehlenn!
%             F = AnisoPassiveTendon(this.stretchfun(0))*pi*this.radfun(0)^2;
%             
%             % Fläche als Funktion des Ortes
%             A = @(y)this.radfun(y).^2*pi;
%             
%             tmrFunc=@(y) (F./A(y)-AnisoPassiveMuscle(this.stretchfun(y)))./...
%                 (AnisoPassiveTendon(this.stretchfun(y))-AnisoPassiveMuscle(this.stretchfun(y)));

            tmrFunc=@(x) exp(-x);
            [ypoints] = getDiscr(this,tmrFunc);
%             ypoints=linspace(0,this.ylen,15);
            plot(ypoints,tmrFunc(ypoints));
            geo= fem.geometry.Belly(ypoints,'Radius',this.radfun,'InnerRadius',0);
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            displ_dir([1 3],geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir([1 3],geo.Elements(13:16,geo.MasterFaces(4,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            velo_dir(2,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            velo_dir(2,geo.Elements(13:16,geo.MasterFaces(4,:))) = true;
%             for k = [1 2 7 8 5 6 11 12]  
%                 pos = geo.Elements(k,geo.MasterFaces(3,:));
%                 velo_dir(1,pos) = true;
%                 %velo_dir_val(1,pos) = 1;
%             end
            velo_dir_val(velo_dir) = 2;
        end
        
        function anull = seta0(~, anull)
            % Direction is y
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        
        function test_ShakerDefaultTMR

            c = ShakerDefault('Stretch','Gauss 0.3');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
            
%             c = ShakerDefault('Stretch','Gauss 0.2');
%             m = c.createModel;
%             c.plotTMR;
%             
%             c = ShakerDefault('Stretch','log');
%             m = c.createModel;
%             c.plotTMR;
%             m.plotGeometrySetup;
%             
% 
%             c = ShakerDefault('Stretch','exp');
%             m = c.createModel;
%             c.plotTMR;
%             m.plotGeometrySetup;
%     
%             c = ShakerDefault('Stretch','linear');
%             m = c.createModel;
%             c.plotTMR;
%             m.plotGeometrySetup;
        end
        
        function test_ShakerDefault
            m   = models.muscle.Model(ShakerDefault);
            m.simulateAndPlot;
%             m = models.muscle.Model(ShakerDefault('Frequency',80));
%             m.simulateAndPlot;
        end
    end
    
end
