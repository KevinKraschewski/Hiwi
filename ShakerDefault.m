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
        maxN;
        maxYLength;
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
            %Default maxYLength to 5
            this.addOption('maxYLength',5);
            %Default N points on yAxis (not the total amount of nodes!)
            this.addOption('maxN',100);
            
            
            % Compute outer shape
            this.ylen = 100;
            k = kernels.GaussKernel(20);
            this.radfun = @(x)10*k.evaluate(this.ylen/2,x)+2;
            
            % Also invokes geo setup
            this.init;
            %Sets the frequency property
            this.frequency=this.Options.Frequency;
%             %Sets the maxN property
%             this.maxN=this.Options.maxN;
%             %Sets the tolerance property
%             this.TOL=this.Options.TOL;
%             %Sets the maxYPoints property
%             this.maxYPoints=this.Options.maxYPoints;
%             %Sets the maxYLength property
%             this.maxYLength=this.Options.maxXLenghth;
            
%             ref = [1.03, 1.4];
%             %Man kann nun ueber 'Gauss Gamma' einen GaußKernel als
%             %Stretchfunktion verwenden, wichtig ist das Leerzeichen
%             %zwischen Gauss und dem gewuenschten Gamma.
%             %Bei groesseren Gamma -> Ungenauigkeiten im maximal Wert, fuehrt zur Unbrauchbarkeit!
%             if length(this.Options.Stretch) >= 5 && strcmp(this.Options.Stretch(1:5),'Gauss')
%                 Gamma=str2double(this.Options.Stretch(6:length(this.Options.Stretch)));
%                 k=kernels.GaussKernel(Gamma);
%                 kexp=kernels.KernelExpansion;
%                 kexp.Kernel=k;
%                 kexp.Centers.xi=1;       % Funktion soll auf [0,1] definiert sein
%                 kexp.Ma=(ref(2)-ref(1)); %Damit insgesamt das Maximum wieder bei ref(2) liegt
%                stretchfun_0_1=@(x) (ref(1)-kexp.evaluate(0))+kexp.evaluate(x); %Der erste Term ist uA da um "Ungenauigkeiten" auszugleichen
%             else
%             switch this.Options.Stretch
%                 case 'exp'
%                     stretchfun_0_1 =@(x) ref(1)+(ref(2)-ref(1)).*x.^4;
%                 case 'log'
%                     stretchfun_0_1 =@(x)((ref(2)-ref(1))/log(51)).*log(50.*(1/50+x))+ref(1);
%                 case 'linear'
%                     stretchfun_0_1 =@(x) (ref(2)-ref(1)).*x+ref(1);
%             end
%             end
%             len = this.ylen;
%             % Transfer the [0,1] argument version into the given domain
%             % length (mirror at half)
%             this.stretchfun = @(y)(y<=len/2) .* (stretchfun_0_1(2*y/len)) ...
%                 + (y>len/2) .* (stretchfun_0_1(1-2*(y-len/2)/len));
            this.maxYLength=this.Options.maxYLength;
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
        
  
%%  Discr due to function values, TMRFunc must be a function handle
     function [ypoints] = getDiscr(this,TMRFunc)
         %Kurze Fehlerbehandlung
         if ~isa(TMRFunc,'function_handle')
             error('The function argument must be of type function_handle');
         end
         if this.maxYLength*this.maxN*2 < this.ylen
            error('Can not reach maxYLength with this amount of points');
         end

         
         %% Vorinitialisierung
         %N=NoElements/4 -> Layer, NoElements etc. werden erst im Anschluss
         %festgelegt
         N=this.maxN;
         N=N/2; %Nur bis zur haelfte
         ypoints=zeros(N,1); 
         distance=zeros(N,1); %Speichere die Abstaende um spaeter an der Mitte zu spiegeln
         j=2; %Erster Punkt ist die 0, fuer den Index in distance und ypoints
         yDefault=linspace(0,this.ylen/2,1000);
         TMRFuncComp=TMRFunc(yDefault); %Vergleichswerte fuer spaeter
         refValue=TMRFunc(0); %Referenzwert um die Toleranz zu checken
         posChange=0; %Wird in der Hauptlogik benoetigt um die variable Laenge der YPunkte zu beruecksichtigen
         
         %% Abschnitt um Genauigkeit in den Vergleichsfunktionswerten zu erreichen
         compare=abs(TMRFuncComp(2:length(TMRFuncComp))-TMRFuncComp(1:length(TMRFuncComp)-1))>this.TOL; %Sollte die Distanz zwischen zwei benachbarten Punkten die Toleranz ueberschreiten -> Refine, damit Rest funktioniert
         compare=find(compare); %Indices der Intervalle mit TOL Ueberschreitung
         while ~isempty(compare) %Solange es Intervalle gibt die nicht passen -> Selbe Logik wie unten
             for k=1 : length(compare)
                 yDefault=[yDefault(1:compare(k)) (yDefault(compare(k))+yDefault(compare(k)+1))/2 yDefault(compare(k)+1:length(yDefault))];
             end
         TMRFuncComp=TMRFunc(yDefault);    
         compare=abs(TMRFuncComp(2:length(TMRFuncComp))-TMRFuncComp(1:length(TMRFuncComp)-1))>this.TOL; %Sollte die Distanz zwischen zwei benachbarten Punkten die Toleranz ueberschreiten -> Refine, damit Rest funktioniert
         compare=find(compare);
         end
         
         
         %% Hauptlogik
         

        for i=1 : length(yDefault);
            
            if(abs(TMRFuncComp(i)-refValue)>this.TOL && i<=length(yDefault)) %Erster Punkt der TOL ueberschreitet -> Vorheriger in ypoints, setzt voraus, dass vorherige Werte diese Toleranz unterschreiten 
                refValue=TMRFuncComp(i-1); % Neuer Vergleichspunkt wird aktueller Funktionswert
                ypoints(j)=yDefault(i-1); %Punkt bei dem TOL noch eingehalten wird, wird hinzugefuegt
                distance(j-1)=ypoints(j)-ypoints(j-1); %Berechne den Abstand benachbarter ypoints Punkte
                j=j+1;
                if j>N 
                %%Alle Punkte verwendet, fertig -> Diverse Abfragen (maxX? gesamtTol?)
                   if(abs(refValue-TMRFunc(this.ylen/2))>this.TOL) %Toleranz insgesamt nicht erreicht 
                        error('Tolerance can not be reached with this amount of points') 
                   end
                   if(abs(max(distance))>this.maxYLength)
                       error('Can not fall below maxYLength with this amount of points');
                   end
                   if (2*abs(ypoints(N)-this.ylen)>this.maxYLength) %distance zwischen letztem Punkt links und erstem rechts
                       error('can not fall below maxYLength with this amount of points');
                   end
                end
               
            elseif i+1==length(yDefault) %Toleranz wird am Ende nicht mehr ueberschritten -> verfuegbarer Punkt wird Mittelpunkt
                ypoints(j)=this.ylen/2;
                ypoints=ypoints(1:j);
                distance(j-1)=ypoints(j)-ypoints(j-1); %Abstaende zwischen den Punkten speichern, fuer spaeter
                distance=distance(1:j-1);
            end
        end
        %%Bis hier sind die Punkte gesetzt um die Toleranz einzuhalten!
        %%weiter sind hier noch Punkte uebrig um zu Verfeinern ->
        %%maxYLength einhalten
        k=find(distance>this.maxYLength);   
        usablePoints=2*N-2*j+1; %Auf dem gesamten Intervall
        while max(distance>this.maxYLength)==1 && length(k)*2 <= usablePoints %Laenge von k sind mind. Anzahl an neuen Punkten auf einer Haelfte
            
            for i=1:length(k) %Iteration ueber die Intervalle in denen maxX ueberschritten wird, zwischen die Punkte einene weiteren setzen.
                ypoints=[ypoints(1:k(i)+posChange); (ypoints(k(i)+posChange)+ypoints(k(i)+1+posChange))/2; ypoints(k(i)+1+posChange:length(ypoints))]; %y veraendert sich nach einem Schleifendurchlauf -> PosChange+1
                distance=[distance(1:k(i)+posChange-1);distance(k(i)+posChange)/2; distance(k(i)+posChange)/2 ; distance(k(i)+1+posChange:length(distance))]; %Abstand entsprechend angepasst -> halbiert dafuer 2mal und posChange+2
                posChange=posChange+1;
                j=j+1; %Fuer das Spiegeln an der Mitte
            end
            usablePoints=usablePoints-2*length(k); %length(k) Punkte auf einer Haelfte gesetzt            
            k=find(distance>this.maxYLength);      
            posChange=0;
        end
        if(max(k)>=1)
            error('Can not fall below maxYLength with this amount of points');
        end

        %Zum "Spiegeln" an der Mitte
        ypointsTemp=zeros(length(ypoints),1);
        distance=flipud(distance);
        ypointsTemp(1)=ypoints(j);
        ypointsTemp(length(ypoints))=ypoints(j)+this.ylen/2;
        for i=2 : length(ypoints)-1
            ypointsTemp(i)=ypointsTemp(i-1)+distance(i-1); %Mit distance wird dafuer gesorgt, dass die Abstaende zwischen den Punkten rechts erhalten bleibt
        end
        
        ypoints=[ypoints(1:j-1) ; ypointsTemp];
        disp(['Maximal distance between yPoints is ',num2str(max(distance)) ,' mm']);
        disp(['Tolerance of ' num2str(this.TOL) ' is reached and ' num2str(length(ypoints)) ' points are used']);
     end
    end
    methods(Access=protected)
        
        function geo = getGeometry(this)
            m  = models.muscle.Model(models.muscle.examples.Belly);
            mu = m.DefaultMu;
            this.TOL=this.Options.TOL;
            this.maxN=this.Options.maxN;
            this.maxYLength=this.Options.maxYLength;
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
            
            S_m = models.muscle.functions.MarkertLawOriginal(mu(5),mu(6));
            AnisoPassiveMuscle = S_m.getFunction;
            S_t = general.functions.CubicToLinear(mu(7),mu(8));
            AnisoPassiveTendon = S_t.getFunction;
            F = AnisoPassiveTendon(this.stretchfun(0))*pi*this.radfun(0)^2;
            
            % Fläche als Funktion des Ortes
            A = @(y)this.radfun(y).^2*pi;
            
            tmrFunc=@(y) (F./A(y)-AnisoPassiveMuscle(this.stretchfun(y)))./...
                (AnisoPassiveTendon(this.stretchfun(y))-AnisoPassiveMuscle(this.stretchfun(y)));
            
            [ypoints] = getDiscr(this,tmrFunc);
            geo = fem.geometry.Belly(ypoints,'Radius',this.radfun,'InnerRadius',0);
            
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            displ_dir([1 3],geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir([1 3],geo.Elements(length(geo.Elements)-3:length(geo.Elements),geo.MasterFaces(4,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            velo_dir(2,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            velo_dir(2,geo.Elements(length(geo.Elements)-3:length(geo.Elements),geo.MasterFaces(4,:))) = true;
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

%             c = ShakerDefault('Stretch','Gauss 0.3');
%             m = c.createModel;
%             c.plotTMR;
%             m.plotGeometrySetup;
            
%             c = ShakerDefault('Stretch','Gauss 0.2');
%             m = c.createModel;
%             c.plotTMR;
%             
            c = ShakerDefault('Stretch','log');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
            

            c = ShakerDefault('Stretch','exp');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
    
            c = ShakerDefault('Stretch','linear');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
        end
        
        function test_ShakerDefault
            m   = models.muscle.Model(ShakerDefault);
            m.simulateAndPlot;
        end
    end
    
end
