%% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!


NoOfSteps = 8; % Wie viele Steps gemacht werden
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

%% Schleife um diverse Optionen zu durchlaufen
for j = 1 : 7
    switch j
        % Nur die Toleranz veraendert sich
        case 1
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2);
            
            % maxYLength veraendert sich
        case 2
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.25,'maxYLength',17.5);
            
            % Toleranz und maxYLength veraendern sich
        case 3
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2);
            
            % Mittelwerte von MuscleTendonRatioGP
        case 4
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2,'MaxYLength',10);
            
            % Mittelwerte, TOL und maxYLength
        case 5
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2);
            
        case 6
            
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2,'constFromTo',[10,40]);
            
        case 7
            c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.2,'constFromTo',[10,40],'maxYLength',17.5);
    end
    
    %% Schleife um NoOfSteps Simulationen durchzufuehren
    for i = 1 : NoOfSteps
        %% Manipulation der Config abhaengig vom Case j
        
        
        
        switch j
            % Nur die Auswirkungen der TOL Veraenderung
            case 1
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',c.TOL-(i-1)*.25);
                
                % Nur die Auswirkungen von maxYLength
            case 2
                              
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.25,'maxYLength',c.maxYLength - (i-1)*2.5);
                % Tol und maxYLength gleichzeitig
            case 3
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',c.TOL-(i-1)*.25,'maxYLength',c.maxYLength - (i-1)*2.5);
                
                
                % Tol und maxYLength gleichzeitig und MW der TMR's
            case 5
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',c.TOL-(i-1)*.25,'maxYLength',c.maxYLength - (i-1)*2.5);
                
                
                % Beides mit Stueckweise Konstanten Muskel/Tendon
            case 6
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',c.TOL-(i-1)*.25,'constFromTo',[10,40],'maxYLength',c.maxYLength - (i-1)*2.5);
                
                % TOL mit Stueckweise Konstant und MW
            case 7
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',c.TOL-(i-1)*.25,'constFromTo',[10,40],'maxYLength',17.5);
                
        end
        
        %% Creating the Model
        m = c.createModel;
        m.T = 150;
        
        
        % Ersetzt die TMR's durch die Mittelwerte der bisherigen
        if j == 4 || j == 5 || j == 7
            % Wichtig, das Attribut wurde in System auf public gesetzt!!
            m.System.MuscleTendonRatioGP(1:size(m.System.MuscleTendonRatioGP,1),:)...
                = repmat(mean(m.System.MuscleTendonRatioGP),...
                size(m.System.MuscleTendonRatioGP,1),1);
            
        end
        
        % Speichert Model und Configuration in eine Cell
        ModelsAndConfig{i,1} = c;
        ModelsAndConfig{i,2} = m;
        
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime(i) = toc(start); % Zeiten auf  ~4.4 GHz
        DF{i} = m.getResidualForces(t,y);
        ElapsedTime{i} = t(length(t)); % Um vorzeitigen Abbruch zu erkennen
        
        %% Saving the Results into a cell
        SimResults{i,1}=t;
        SimResults{i,2}=y;
        
        % Nur ein Durchgang, falls nur die MW betrachtet werden sollen
        if j == 4
            break
        end
        
        
        
    end
    
    %% Speichert die meisten Variablen in eine Datei "Results".
    save(['Results' num2str(j) '.mat'],'SimResults','SimulationTime','ElapsedTime','NoOfSteps','DF','ModelsAndConfig');
end