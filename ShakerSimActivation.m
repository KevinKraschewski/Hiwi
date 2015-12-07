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
SimResults_Time = cell(NoOfSteps,2);

% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(NoOfSteps,1);

% Speichert die untersch. Kraefte auch in einer Cell
DF = cell(NoOfSteps,1);

%% Schleifen fuer den Aktivierungs/Frequenzen Part
for i = 1 : 8
    
    
    parfor j = 1 : 20
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(50-(i-1)*5,12/10,1-(j-1)*(0.05),30);
        m.dt = .25
        m.DefaultMu(2) = 10;
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime = toc(start); % Zeiten auf  ~4.4 GHz
        DF{j} = m.getResidualForces(t,y);
        ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        %% Saving the Results into a cell
        SimResults_Time{j}=t;
        SimResults_Y{j}=y;
        % Speichert Model und Configuration in einer Cell
        Config{j} = c;
        Model{j} = m;
    end
    %% Speichert die meisten Variablen in eine Datei "Results".
    save(['SmallActivationChange' num2str(i+3) '_12mm' '.mat'],'SimResults_Time','SimResults_Y','SimulationTime','ElapsedTime','DF','Config','Model');
end

%% Schleife fuer den Amplituden/Frequenz Part
for i = 1 : 8
    parfor j = 1 : 18
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(50-(i-1)*5,(12-(j-1)/2)/10,1,30);
        m.dt = .25
        m.DefaultMu(2) = 10;
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime = toc(start); % Zeiten auf  ~4.4 GHz
        DF{j} = m.getResidualForces(t,y);
        ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        %% Saving the Results into a cell
        SimResults_Time{j}=t;
        SimResults_Y{j}=y;
        % Speichert Model und Configuration in einer Cell
        Config{j} = c;
        Model{j} = m;
    end
    
    save(['SmallAmplitudeChange' num2str(i) '.mat'],'SimResults_Time','SimResults_Y','SimulationTime','ElapsedTime','DF','Config','Model');
    
end