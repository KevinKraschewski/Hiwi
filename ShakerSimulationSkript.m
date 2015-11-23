%% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!


NoOfSteps = 1; % Wie viele Steps gemacht werden
TolStep = .025; % Veraenderung der Toleranz mit jedem Schritt

ModelsAndConfig = cell(NoOfSteps,2); % Spalte 1 -> Configs, Spalte 2 -> Models

% Speichert im Eintrag i des Vektors die Dauer der Simulation fuer die i.
% Iteration
SimulationTime = zeros(NoOfSteps,1);

% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults = cell(NoOfSteps,2);

% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(NoOfSteps,1);

% Speichert die untersch. Kraefte auch in einer Cell
DF = cell(NoOfSteps,1);

% Anfangs relative Toleranz fuer den ODE Solver
RelTol = .1;

%% Schleife um diverse Optionen zu durchlaufen
for j = 1 : 3
    c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.25,...
        'constFromTo',[10,40],'maxYLength',20);
    tolerance = .25;
    laenge = 20;
    
    
    %% Schleife um NoOfSteps Simulationen durchzufuehren
    for i = 1 : NoOfSteps
        
        
        % Veraendert die Feinheit der Zerlegung und die relative Toleranz
        % des ODE Solvers in jedem Schritt
        if j == 1
            
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',tolerance-(i-1)*TolStep,...
                'constFromTo',[10,40],'maxYLength',laenge - (i-1)*2.5);
            
            m = c.createModel;
            m.T = 150;
            RelTol = RelTol/10;
            m.ODESolver.RelTol = RelTol;
        end
        
        % Bessere Gauss Integration
        if j == 2
            
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',tolerance-(i-1)*TolStep,...
                'constFromTo',[10,40],'maxYLength',laenge - (i-1)*2.5);
            m = c.createModel;
            m.T = 150;
            m.setGaussIntegrationRule(5);
            
        end
        
        % Auswirkung von beidem
        if j == 3
            
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',tolerance-(i-1)*TolStep,...
                'constFromTo',[10,40],'maxYLength',laenge - (i-1)*2.5);
            m = c.createModel;
            m.T = 150;
            m.setGaussIntegrationRule(5);  
            RelTol = RelTol/10;
            m.ODESolver.RelTol = RelTol;            
            
        end
        % Ersetzt die TMR's durch die Mittelwerte der bisherigen
        % Wichtig, das Attribut wurde in System auf public gesetzt!!
        m.System.MuscleTendonRatioGP(1:size(m.System.MuscleTendonRatioGP,1),:)...
            = repmat(mean(m.System.MuscleTendonRatioGP),...
            size(m.System.MuscleTendonRatioGP,1),1);
        
        % Kleine Anpassung um an Stellen wirklich nur Tendon/Muscle zu
        % haben
        m.System.MuscleTendonRatioGP = (m.System.MuscleTendonRatioGP - ...
            min(m.System.MuscleTendonRatioGP(:)))/(max(m.System.MuscleTendonRatioGP(:)...
            -min(m.System.MuscleTendonRatioGP(:))))./(10^(8));
        
        
        % Speichert Model und Configuration in einer Cell
        ModelsAndConfig{i,1} = c;
        ModelsAndConfig{i,2} = m;
        
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime(i) = toc(start); % Zeiten auf  ~4.4 GHz
        DF{i} = m.getResidualForces(t,y);
        ElapsedTime{i} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        %% Saving the Results into a cell
        SimResults{i,1}=t;
        SimResults{i,2}=y;
        
    end
    
    %% Speichert die meisten Variablen in eine Datei "Results".
    save(['SmallTMR' num2str(j) '.mat'],'SimResults','SimulationTime','ElapsedTime','NoOfSteps','DF','ModelsAndConfig');
end