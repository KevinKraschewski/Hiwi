classdef GetRadiusMuscle < models.muscle.AMuscleConfig

% For given external Forces/Pressure calculates a fitting radius of the
% Muscle to given inner Radius r
%properties
%    TendonRadius;
%    Pressure;
%    Xpkts;
%end

%
%%
% Für jeden Punkt in xgrid:
% Hole die Spannung  an dem Punkt für den Muskel und für die Sehne (eval
% Funktion in Dynamics  (Daniel))
%

    %Hole die Spannung an dem Punkt fuer den Muskel und fuer die Sehne, am
    %Besten Vektoren!
    
    
% Ueber GG der Kräfte holt man sich ueber
% SQRT(Stress_Tendon*A_Tendon/Stress_Muscle*PI)=R_Muscle
% 
methods(Static)
        % Soll zwei Vektoren mit der Spannung von Sehne/Muskel zurueck
        % geben
        function [S_t,S_m] = getStress(this,xgrid) %Evtl. mehr Argumente -> Eval
        sys = this.fsys;
        mc = sys.Model.Config;
        fe_pos = mc.FEM;
        geo = fe_pos.Geometry;
        fe_press = mc.PressureFEM;
            num_u_glob = geo.NumNodes*3;
            idx_u_elems_local = sys.idx_u_elems_local;
            num_elements=numel(xgrid);
            for m = 1:num_elements
        % 1:num_u_glob is all u
        elemidx_u = idx_u_elems_local(:,:,m); 
        % num_u_glob next ones are all v
        elemidx_v = num_u_glob + elemidx_u;
        
        u = uvwcomplete(elemidx_u);
        w = uvwcomplete(idx_p_elems_global(:,m));
        
        if havefibretypes 
            ftwelem = fibretypeweights(:,:,m)*FibreForces;
        end

        integrand_u = zeros(3,dofsperelem_u);
        %Iteration ueber die Pkte in einem Element
        for gp = 1:num_gp 

            % Evaluate the pressure at gauss points
            p = w' * fe_press.Ngp(:,gp,m);

            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);

            if any(isnan(u(:)))
                fprintf('NaNs in models.muscle.Dynamics#evaluateCoreFun at t=%g! Have a look.\n',t);
                keyboard;
            end
            % Deformation gradient
            F = u * dtn;
            C = F'*F;
           
            %% Isotropic part (Invariant I1 related)
%             I1 = sum(sum((u'*u) .* (dtn*dtn')));
            I1 = C(1,1) + C(2,2) + C(3,3);
            
            %% Compile tensor
            P = mooneyrivlin_ic_const(gp,m)*Id3 + p*inv(F)' + 2*(c10(gp,m) + I1*c01(gp,m))*F ...
                - 2*c01(gp,m)*F*C;
            
            %% Anisotropic part (Invariant I4 related)
            if havefibres
                fibrenr = (m-1)*num_gp + gp;
                fibres = sys.a0Base(:,:,fibrenr);
                lambdaf = norm(F*fibres(:,1));
                
                % Get weights for tendon/muscle part at current gauss point
                if hastendons
                    tendonpart = tmrgp(gp,m);
                    musclepart = 1-tendonpart;
                end
                if havefibretypes
                    alpha = musclepart*ftwelem(gp);
                else
                    alpha = musclepart*alphaconst;
                    %[t sys.MuscleTendonRatiosGP(gp,m) alpha]
                end
                passive_aniso_stress = 0;
                % Using > 1 is deadly. All lambdas are equal to one at t=0
                % (reference config, analytical), but numerically this is
                % dependent on how precise F and hence lambda is computed.
                % It is very very close to one, but sometimes 1e-7 smaller
                % or bigger.. and that makes all the difference!
                if lambdaf > 1.0001
                    passive_aniso_stress = musclepart*anisomusclefun(lambdaf);
                    if hastendons
                        passive_aniso_stress = passive_aniso_stress + tendonpart*anisotendonfun(lambdaf);
                    end
                end
                
                % Using a class-subfunction is 20% slower!
                % So: function handle
                fl = flfun(lambdaf);
                gval = passive_aniso_stress + (Pmax/lambdaf)*fl*alpha;
                P = P + gval*F*sys.a0oa0(:,:,fibrenr);
                
                %% Cross-fibre stiffness part
                if usecrossfibres
                    lambdaf = norm(F*fibres(:,2));
                    if lambdaf > .999
                        g1 = (b1cf/lambdaf^2)*(lambdaf^d1cf-1);
                        P = P + g1*F*sys.a0oa0n1(:,:,fibrenr);
                    end
                    lambdaf = norm(F*fibres(:,3));
                    if lambdaf > .999
                        g2 = (b1cf/lambdaf^2)*(lambdaf^d1cf-1);
                        P = P + g2*F*sys.a0oa0n2(:,:,fibrenr);
                    end
                end
                
                %% Check if change rate of lambda at a certain gauss point should be tracked
                % (corresponds to a spindle location in fullmodels.muscle.Model)
                if ~isempty(ldotpos)
                    k = find(ldotpos(1,:) == m & ldotpos(2,:) == gp);
                    if ~isempty(k)
                        Fdot = uvwcomplete(elemidx_v) * dtn;
                        this.lambda_dot(k) = (F*fibres(:,1))'*(Fdot*fibres(:,1))/lambdaf;
                    end
                end
            end
            
           %%  Viscosity - currently modeled as extra A linear system
%             if visc > 0
%                 v = uvw_full(elemidx_v);
%                 P = P + visc * v * dtn;
%             end

            %% Assembly part I - sum up contributions from gauss points
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m,gp);

            integrand_u = integrand_u + weight * P * dtn';
        end % end of gauss point loop
            end
        end
        
        %Gives back a radius vector for the muscle for given stresses and
        %tendon radius, 
        %"xgrid" just needed for the loop indices -> now just Vector ops,
        %no loop
        %Test mit S_t,S_m,r_t Zeilenvektor
function [r_m] = getRadius(S_t,S_m,r_t)
  % Make sure its a row vector
        S_t = reshape(S_t,1,[]);
        a   = length(S_t);
        S_m = reshape(S_m,1,[]);
        b   = length(S_m);
        r_t = reshape(r_t,1,[]);
        
        %Fuer S_t gleiche Laenge S_m
        
        if (length(r_t) ~= length(S_m))
            r_t(length(r_t):length(S_m)) = r_t(length(r_t)); %sinnvoll?
        end
        
        %%Eingaben untersch längen, sollte bei getStress nicht vorkommen,
        %%da jeweils Stress fuer jedes Element aus x zurueck gibt.
        
        %Nicht zufrieden!!!
        %h=[length(S_t),length(S_m),length(r_t)];
        
        %max_Laenge=max(h);
        
        %Es treten nicht alle Faelle gleichzeitig auf, eine Abfrage
        %einsparen?
        %if(max_Laenge > length(S_m))
        %    S_m = interp1(1:length(S_m), S_m, linspace(1, length(S_m), max_Laenge))
        %end
        %if(max_Laenge > length(S_t))
        %    S_t = interp1(1:length(S_t), S_t, linspace(1, length(S_t), max_Laenge))
        %end
        %if(max_Laenge > length(r_t))
        %    r_t = interp1(1:length(r_t), r_t, linspace(1, length(r_t), max_Laenge))
        %end
        
        %% Calculate r_m
        
 %      A_t = zeros(numel(xgrid),1);
 %      r_m = zeros(numel(xgrid),1);
 %           for i=1:numel(xgrid)
 %               A_t(i)=PI*r_t(i)^2; %Kreisfoermige Flaeche
 %               r_m(i)=sqrt(S_t(i)*A_t(i)/(S_m(i)*PI)); %Ueber Vektor Operationen?!
 %           end
 
 
 
 
 
            A_t=pi.*r_t.^2; %Kreisfoermige Flaeche
            r_m=sqrt(S_t.*A_t./(S_m.*pi));
       end

% Konstruiere die Geometrie, mittels vorhandenen Methoden
% 
function geo = construct(xgrid,r_t)
     if isscalar(r_t)
         r_t=ones(numel(xgrid),1)*r_t; %Spaltenvektor
     end
%    [S_t,S_m]   = getStress(xgrid); %Noch nicht existent
     S_t=linspace(1,max(xgrid),length(xgrid)).*5; %Wird nicht mehr benoetigt, wenn getStress funktioniert
     S_m=linspace(1,max(xgrid),length(xgrid)).*5^(-2);
     r_m = GetRadiusMuscle.getRadius(S_t,S_m,r_t); %Vector mit gl Dim wie xgrid, Zeilenvektor
    
     geo = fem.geometry.Belly(xgrid,'Radius',r_m'-r_t,'InnerRadius',r_t,'Gamma',2);
    
end
    function GetRadiusTest
     geo =   GetRadiusMuscle.construct([0 .5 .6 .7 .8 .9 1 1.5 2 2.4 2.5 2.7 3 3.5 4 5 6 6.2 7.1 8 9],.5);
     
     geo.plot %Plot falsch? bekomme auf jedem Pkt ein r_m, wächst in der Konstellation zu Beginn schneller
     % und am Ende langsamer.. plot liefert aber nur an ein paar Punkten
     % die "richtige" Konstruktion?! --> needs to be fixed. (Funktioniert,
     % wenn Vektoren nicht zu groß werden?!)
    end
end
end

% Radius muscle is Outer Radius-Inner Radius
%
%

%%methods(Static)           
%            function Geom = Funktion(InnerR, np)            
%            Geom = fem.geometry.Belly(np,'InnerRadius',InnerR,'Gamma',2);
%            
%            plot(Geom)
%            end
%end
    

%% Teile aus evaluate@Dynamics
%       sys = this.fsys;
%       mc = sys.Model.Config;
%       fe_pos = mc.FEM;
%       fe_press = mc.PressureFEM;
% m ist der Index, der aeusseren Schleife (NumElements von geo dh.
% geo.NumElements)
%       w = uvwcomplete(idx_p_elems_global(:,m));

%       Nur ein Pkt xpkt -> gp=xpkt
%       p = w' * fe_press.Ngp(:,xpkt,m);
