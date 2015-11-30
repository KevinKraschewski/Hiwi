%% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!


NoOfSteps = 1; % Wie viele Steps gemacht werden

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

% Schleife zum aendern der maximal Aktivierung und der
% Vollaktivierungsdauer
for j = 1 : 4
    
    
    % Schleife ueber die unterschiedlichen Frequenz/Amplituden
    % Kombinationen
    for i = 1 : 4
        
        switch i
            case 1
                
                [c,m] = ShakerDefaultFuerSim.createTestingConfig(15,4/10,1/j^2,10*pi);
                
            case 2
                
                [c,m] = ShakerDefaultFuerSim.createTestingConfig(15,12/10,1/j^2,10*pi);
                
            case 3
                
                [c,m] = ShakerDefaultFuerSim.createTestingConfig(50,4/10,1/j^2,10*pi);
                
            case 4
                
                [c,m] = ShakerDefaultFuerSim.createTestingConfig(50,12/10,1/j^2,10*pi);
                
        end
        
        m.DefaultMu(2) = 5*j;
        
        % Speichert Model und Configuration in einer Cell
        ModelsAndConfig{1,1} = c;
        ModelsAndConfig{1,2} = m;
        
        
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime = toc(start); % Zeiten auf  ~4.4 GHz
        DF{1} = m.getResidualForces(t,y);
        ElapsedTime{1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        %% Saving the Results into a cell
        SimResults{1,1}=t;
        SimResults{1,2}=y;
        
        %% Speichert die meisten Variablen in eine Datei "Results".
        save(['Activation' num2str(i) '_' num2str(j) '.mat'],'SimResults','SimulationTime','ElapsedTime','DF','ModelsAndConfig');
    end
end