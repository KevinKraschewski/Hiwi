classdef TMR<models.muscle.AMuscleConfig
    methods(Static)
        %TMR.buildGeo erzeugt Geometrie mit den Parameter:
        %len=laenge des Muskel/Sehnen Konstrukts
        %r_t = Tendonradius
        %r_m = Muskelradius
        function geo = buildGeo(len,r_t,r_m)
            geo = fem.geometry.Belly(len,'InnerRadius',r_t,'Radius',r_m);
        end
        
        %%TMR.getDistortionFunction(y,trend) gibt eine Funktion lambda zurueck, die einen Verzerrungsverlauf zwischen y(0) und
        %%y(1) interpoliert (x aus [0,1]!!)
        %mit 'trend' laesst sich der Verlauf als exponentialfunktion
        %('exp'), als logarithmus ('log') oder als linear ('linear') einstellen.
        
        
        %Plus r_a(x)?
        function [lambda] = getDistortionFunction(y,trend)
            if ~isstring(trend)
                error('Please use a string for the second input argument')
            end
            if ~strcmp(trend,'exp') || ~strcmp(trend,'log') || ~strcmp(trend,'linear')
                    error('trend must be either "log", "linear" or "exp"')
            end
            
            switch trend
                case 'exp'
                    lambda =@(x) y(1).*exp(x.*(log(y(2)/y(1))));
                case 'log'
                    
                    lambda =@(x) ((y(2)-y(1))/log(2)).*log(1+x)+y(1); 
                case 'linear'
                    lambda =@(x) (y(2)-y(1)).*x+y(1);
            end    
        
        end
        %%Gibt eine Stressfunktion zurueck.
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
    
    end
end