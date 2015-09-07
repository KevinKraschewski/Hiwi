classdef GetRadiusMuscle_structured < models.muscle.AMuscleConfig

% For given external Forces/Pressure calculates a fitting radius of the
% Muscle to given inner Radius r


%
%%
% Für jeden Punkt in xgrid:
% Hole die Spannung  an dem Punkt für den Muskel und für die Sehne (eval
% Funktion in Dynamics  (Daniel))
%

    %Hole die Spannung an dem Punkt fuer den Muskel und fuer die Sehne, am
    %Besten Vektoren!
    
    
% Ueber GG der Kräfte holt man sich ueber
% 
methods(Static)
        % Soll zwei Vektoren mit der Spannung von Sehne/Muskel zurueck
        % geben
        function [S_t,S_m] = getStress(this,xgrid) %Evtl. mehr Argumente -> Eval

        end

        %%GetRadius Methode: Liefert fuer gegebenen Stress und
        %%Tendonradiusverlauf einen passenden Muskelradiusvektor zurueck
        function [r_m] = getRadius(S_t,S_m,r_t)
        % Make sure its a row vector
        S_t = reshape(S_t,1,[]);
        S_m = reshape(S_m,1,[]);
        r_t = reshape(r_t,1,[]);
        
        %Fuer S_t gleiche Laenge S_m, sollte ueber kommende getStress
        %Methode gewaehrleistet sein
        
        %% Calculate r_m
        

           r_m = zeros(numel(S_t),1);
           A_t=pi.*r_t.^2; %Kreisfoermige Flaeche
           
            for i=2:numel(S_t) %i=2, da fuer i=1 nur Tendon , fest!
              
                A_m = r_m(i-1)^2*pi; %Flaeche Muskel von i-1
                f_m = S_m(i-1)*A_m; %Kraft Muskel von i-1
                
                
                %Berechnung ueber GG zwischen Punkt i und i-1
                r_m(i)=sqrt(abs((S_t(i-1)*A_t(i-1)-S_t(i)*A_t(i)+f_m))/(S_m(i)*pi));
            end
 
       end

% Konstruiere die Geometrie, mittels vorhandenen Methoden
% 
    function geo = construct(xgrid,r_t)
        if isscalar(r_t)
             r_t=ones(numel(xgrid),1)*r_t; %Spaltenvektor
        end
     
%    [S_t,S_m]   = getStress(xgrid); %Noch nicht existent
     S_t=linspace(3,max(xgrid),length(xgrid)).*100; %Wird nicht mehr benoetigt, wenn getStress funktioniert
     S_m=linspace(1,max(xgrid),length(xgrid)); %Dito
     
        if (length(r_t) ~= length(S_m))
           r_t(length(r_t):length(S_m)) = r_t(length(r_t)); %sinnvoll den Vekrot einfach weiterzufuehren?
        end
           
        r_m = GetRadiusMuscle_structured.getRadius(S_t,S_m,r_t); %Vector mit gl Dim wie xgrid, Zeilenvektor    
        geo = fem.geometry.Belly(xgrid,'Radius',r_m'-r_t,'InnerRadius',r_t,'Gamma',3);
     
    end
    


    function GetRadiusTest
        xgrid  =[0 .5 1 1.5 2 2.4 2.5 2.7 3 3.5 4 5 6 6.2 7.1 8 9 9.3 9.5 10 10.9 12 12.5 13 13.4 14 14.3 14.6 14.7 14.9 15 16 17];
        r_t    =linspace(1,.01,10);
 %       r_t    = 1;
        geo    =   GetRadiusMuscle_structured.construct(xgrid,r_t);
        geo.plot %Plot "schlecht" fuer grosse Werte in xgrid (Falsche Berechnung oder Plot?)
    end
end


end