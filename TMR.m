classdef TMR
    methods(Static)
        %%TMR.TestTMR testet die bisherigen Funktionen und plottet die
        %%Auswertungen auf dem Intervall [0,1].
        function TestTMR
            %Einige Parameter, fuer Erklaerungen siehe entsprechende
            %Funktionen
            m  = models.muscle.Model(models.muscle.examples.Belly);
            mu = m.DefaultMu;
            x=linspace(0,1,250);
            y=[1.03,1.4];
            len=1;          
            r_t=@(x) 1;
            r_m=@(x) 2;
            constUntil=0.1;
            
            
            %Distortionfunctions
            figure
            hold on            
            Distortion1=TMR.getDistortionFunction(y,'exp');
            Distortion2=TMR.getDistortionFunction(y,'log');
            Distortion3=TMR.getDistortionFunction(y,'linear');
            plot(x,Distortion1(x))
            plot(x,Distortion2(x));
            plot(x,Distortion3(x));
            legend('exp','log','linear','Location','northwest');
            title('Different Distortionfunction interpolations');
            xlabel('position');
            ylabel('Strain');
            hold off
            
             
             %TendonMuscleRatio plots fuer verschiedene Interpolationenen
             %der Verzerrung
             % 1 = Nur Tendon
             % 0 = Nur Muscle
             figure
             hold on
             title('TMR');
             TenMR1=TMR.getTendMuscRatio(mu,r_t,r_m,y,'exp',constUntil);
             TenMR2=TMR.getTendMuscRatio(mu,r_t,r_m,y,'log',constUntil);
             TenMR3=TMR.getTendMuscRatio(mu,r_t,r_m,y,'linear',constUntil);
             plot(x,TenMR1(x));
             plot(x,TenMR2(x));
             plot(x,TenMR3(x));
             ylabel('Ratio');
             xlabel('Position');
             legend('Dist=exp','Dist=log','Dist=linear','Location','northeast')
             hold off
             
%            Stressfunctions
             Stress=TMR.getStressFunction;

             figure
                
             [X,Y]=meshgrid(Distortion1(x),TenMR1(x));
             mesh(X,Y,Stress(X,Y));
             zlabel('Stress');
             ylabel('Ratio');
             xlabel('Strain');
             title('Stress with exp strain');

             
             
             figure
             
             
             [X,Y]=meshgrid(Distortion2(x),TenMR2(x));
             mesh(X,Y,Stress(X,Y));
             zlabel('Stress');
             ylabel('Ratio');
             xlabel('Strain');
             title('Stress with log strain');
             
             figure
             
             [X,Y]=meshgrid(Distortion3(x),TenMR3(x));
             mesh(X,Y,Stress(X,Y));
              xlabel('Strain');
              ylabel('Ratio');
              zlabel('Stress');
              title('Stress with linear strain');
             hold off

             
             figure
             hold on
             a=@(x) r_t(x).*TenMR1(x);
             plot(x,a(x));
             b=@(x) r_m(x).*(1-TenMR1(x));
             plot(x,b(x));
             legend('Radius Tendon','Radius Muscle');
             title('Development of the Radii');
             xlabel('position')
             ylabel('radius [mm]');
             hold off
             
             Geom=TMR.buildGeo(len,b(x),a(x));
             Geom.plot;        
            
        end
        
        
        
        
        %TMR.buildGeo erzeugt Geometrie mit den Parameter:
        %len=laenge des Muskel/Sehnen Konstrukts
        %r_t = Tendonradius
        %r_m = Muskelradius
        function geo = buildGeo(len,r_t,r_m)
            geo = fem.geometry.Belly(len,'Radius',r_m,'InnerRadius',r_t);
        end
        
        %%TMR.getDistortionFunction(y,trend) gibt eine Funktion lambda zurueck, die einen Verzerrungsverlauf zwischen y(0) und
        %%y(1) interpoliert und diese plottet.
        %!!x aus [0,1]!!
        %mit 'trend' laesst sich der Verlauf als exponentialfunktion
        %('exp'), als logarithmus ('log') oder als linear ('linear') einstellen.
        
        
        function [lambda] = getDistortionFunction(y,trend)
            if ~ischar(trend)
                error('Please use a string for the second input argument')
            end
            if ~strcmp(trend,'exp') && ~strcmp(trend,'log') && ~strcmp(trend,'linear')
                    error('trend must either be "log", "linear" or "exp"')
            end
            
            switch trend
                case 'exp'
                    lambda =@(x) y(1)+(y(2)-y(1)).*x.^4;
                case 'log'                    
                    lambda =@(x)((y(2)-y(1))/log(51)).*log(50.*(1/50+x))+y(1);
                case 'linear'
                    lambda =@(x) (y(2)-y(1)).*x+y(1);
            end    

            
        end
        
        %%TMR.getStressFunction(mu) gibt eine Stressfunktion zurueck.
        %Eingabeparameter ist der Paramter mu
        %Kein inputArgument -> DefaultMu
        function[Stress]=getStressFunction(mu)
            if nargin < 1
                m  = models.muscle.Model(models.muscle.examples.Belly);
                mu = m.DefaultMu;
            end
            
        S_m = models.muscle.functions.MarkertLawOriginal(mu(5),mu(6));
        S_m = S_m.getFunction;
        S_t = general.functions.CubicToLinear(mu(7),mu(8));
        S_t = S_t.getFunction;
        
        Stress=@(lambda,ratio) S_m(lambda).*(1-ratio)+S_t(lambda).*ratio;
            
        end
        
        %%TMR.getTendMuscRatio(Stress,r_t,r_m) soll den
        %%Sehnen-Muskelquotienten fuer eine gegebene Stressfunktion, sowie
        %%inneren und aeusseren Radius der Geometrie zurueckgeben.
        %Paramter:
        %mu ist der Parameter der Funktion getStressFunction.
        %r_t,r_m sind Funktionen abhaengig vom Ort x.
        %y gibt das Verzerrungsintervall vor (wie bei
        %getDistortionFunction).
        %trend liefert die Interpolationsfunktion fuer
        %getDistortionFunction.
        
        %constUntil ist der Parameter, welcher bestimmt bis zu welchem
        %Punkt es nur Sehne geben soll
        function [TendonMuscleRatio]=getTendMuscRatio(mu,r_t,r_m,y,trend,constUntil)

            Stress=TMR.getStressFunction(mu);
            Distortion=TMR.getDistortionFunction(y,trend);
            
            %Berechne zunaechst die Kraft die an der Sehne aufgebracht
            %werden muss, um die von y(1) vorgegebene Verzerrung zu erreichen
            F=Stress(Distortion(0),1)*pi*r_t(0)^2;
            
            %Flaeche als Funktion des Ortes
            A=@(x) (r_t(x)+r_m(x)).^2*pi;
            
            %Der Quotient ausgehend von F/A=Stress
            
            
            TendonMuscleRatio=@(y)(y>=constUntil).*(F-Stress(Distortion(y-constUntil),0))./(A(y).*(Stress(Distortion(y-constUntil),1)-Stress(Distortion(y-constUntil),0)))...
                +(y<constUntil);
               
            
        end
    
    end
end
