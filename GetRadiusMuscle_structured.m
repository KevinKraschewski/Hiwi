%% Mittels "GetRadiusMuscle_Structured.GetRadiusTest" laesst sich die Funktion der Klasse testen.
%Dabei wird ein Vektor xgrid erzeugt, der die einzelnen Punkte in eine
%Raumdimension repraesentiert. Hinzu wird ein Radiusverlauf, r_t, vorgegeben und
%die Methode "construct(xgrid,r_t)" erzeugt mit diesen Werten eine passende
%Geometrie mittels Kraeftegleichgewichten an den einzelnen Punkten in
%xgrid.


classdef GetRadiusMuscle_structured < models.muscle.AMuscleConfig
    

methods(Static)
        %%Statische Methode, die die construct Methode ausfuehrt.
        %Methode gibt Parameter vor.
        function GetRadiusTest
               
                xgrid  = linspace(0,20,25);
                r_t    = 5;
                geo    = GetRadiusMuscle_structured.construct(xgrid,r_t);
                geo.plot
        end
        
%%GetRadiusMuscle_structured.construct(xgrid,r_t,lambda) erzeugt die Geometrie
%Die Methode verwendet dabei die Methoden "getRadius" und "Belly"
%xgrid liefert den Intervallvektor
%r_t ist der Verlauf des Tendonradius
%lambda ist die sich einstellende Verzerrung (der Verlauf)

    function geo = construct(xgrid,r_t,lambda)
        if nargin < 3
            lambda = linspace(1.04, 1.3,numel(xgrid));
        end
        
        if isscalar(r_t)
             r_t=ones(numel(xgrid),1)'*r_t;
        end
     
        [S_t,S_m]=GetRadiusMuscle_structured.getStress(lambda);
        
     
        if (length(r_t) ~= length(S_m))
           r_t(length(r_t):length(S_m)) = r_t(length(r_t)); 
        end
           
        r_m = GetRadiusMuscle_structured.getRadius(S_t,S_m,r_t);
        geo = fem.geometry.Belly(xgrid,'Radius',r_m-r_t,'InnerRadius',r_t,'Gamma',2);
     
    end
        
        %%getStress Methode, holt fuer Verzerrung Lambda und Parameter mu
        %%die Stressverlaeufe fuer tendon und muscle
        
        function [S_t,S_m] = getStress(lambda,mu)
           
            
        if nargin < 2
           m  = models.muscle.Model(models.muscle.examples.Belly);
           mu = m.DefaultMu;
             if nargin < 1
                    error('Please use at least the parameter for the distortion.');
             end
        end   
            
        S_m = models.muscle.functions.MarkertLawOriginal(mu(5),mu(6));
        S_m = S_m.getFunction;
        S_m = S_m(lambda);
        S_t = general.functions.CubicToLinear(mu(7),mu(8));
        S_t = S_t.getFunction;
        S_t = S_t(lambda);
        end

        %%GetRadius Methode: Liefert fuer gegebenen Stress und
        %%Tendonradiusverlauf einen passenden Muskelradiusvektor zurueck
        %S_t ist der Stressverlauf von tendon
        %S_m ist der Stressverlauf von muscle
        %r_t ist der Radiusverlauf von tendon
        
        function [r_m] = getRadius(S_t,S_m,r_t)
        %% Make sure its a row vector
        S_t = reshape(S_t,1,[]);
        S_m = reshape(S_m,1,[]);
        r_t = reshape(r_t,1,[]);
        
        
        
        %% Calculate r_m
           A_t=pi.*r_t.^2; %Kreisfoermige Flaeche          
           %Fliegt um die Ohren!
           r_m = sqrt(S_t.*A_t./(S_m.*pi));           
 
        end

    
end


end