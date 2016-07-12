%% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!


NoOfSteps = 20; % Wie viele Steps gemacht werden

% Speichert im Eintrag i des Vektors die Dauer der Simulation fuer die i.
% Iteration
SimulationTime = cell(NoOfSteps,1);

% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(NoOfSteps,1);

Config = cell(NoOfSteps,1);

Model = cell(NoOfSteps,1);

SimResults_Y = cell(NoOfSteps,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(NoOfSteps,1);

% Speichert die untersch. Kraefte auch in einer Cell
DF = cell(NoOfSteps,1);



%% Schleifen fuer den Aktivierungs/Frequenzen Part
for j = 1 : 13
    
    
    parfor i = 1 : 20
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(90-(j-1)*2.5,12/10,1-(i-1)*(0.05),1);
        m.dt = .25;
        m.DefaultMu(2) = 5;
        %% The Simulation Part with saving the used time into a vector
        start = tic;
        [t,y] = m.simulate;
        SimulationTime{i,1} = toc(start); % Zeiten auf  ~4.4 GHz
        %     DF{i,1} = m.getResidualForces(t,y);
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        if t(end) < 175
            m.ODESolver.RelTol = 10^(-6);
            start = tic;
            [t,y] = m.simulate;
            SimulationTime{i,1} = toc(start); % Zeiten auf  ~4.4 GHz
            %         DF{i,1} = m.getResidualForces(t,y);
            ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        end
        
        if t(end) < 175
            m.ODESolver.RelTol = 10^(-10);
            start = tic;
            [t,y] = m.simulate;
            SimulationTime{i,1} = toc(start); % Zeiten auf  ~4.4 GHz
            %         DF{i,1} = m.getResidualForces(t,y);
            ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        end
        %% Saving the Results into a cell
        SimResults_Time{i,1}=t;
        SimResults_Y{i,1}=y;
        % Speichert Model und Configuration in einer Cell
        Config{i,1} = c;
        Model{i,1} = m;
        i
    end
    j
    %% Speichert die meisten Variablen in eine Datei "Results".
    save(['D:\Hiwi_SimErgebnisse\SmallActivationChange90to60' num2str(j) '_12mm.mat'],'SimResults_Time','SimResults_Y','SimulationTime','ElapsedTime','DF','Config','Model');
end
% save('D:\Hiwi_SimErgebnisse\SmallActivationChange 5 to 38 no activation_12mm.mat','SimResults_Time','SimResults_Y','SimulationTime','ElapsedTime','Config','Model');
% %% Schleife fuer den Amplituden/Frequenz Part
% for i = 1 : 8
%     parfor j = 1 : 18
%         [c,m] = ShakerDefaultFuerSim.createTestingConfig(50-(i-1)*5,(12-(j-1)/2)/10,1,30);
%         m.dt = .25
%         m.DefaultMu(2) = 10;
%         %% The Simulation Part with saving the used time into a vector
%         start = tic;
%         [t,y] = m.simulate;
%         SimulationTime = toc(start); % Zeiten auf  ~4.4 GHz
%         DF{j} = m.getResidualForces(t,y);
%         ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%
%         %% Saving the Results into a cell
%         SimResults_Time{j}=t;
%         SimResults_Y{j}=y;
%         % Speichert Model und Configuration in einer Cell
%         Config{j} = c;
%         Model{j} = m;
%     end
%
%     save(['SmallAmplitudeChange' num2str(i) '.mat'],'SimResults_Time','SimResults_Y','SimulationTime','ElapsedTime','DF','Config','Model');
%
% end