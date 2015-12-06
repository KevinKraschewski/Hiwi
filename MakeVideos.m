% close all
% frequency = zeros(70*20,1);
% maxAmp = zeros(20,1);
% for i = 1 : 70
%     load(['SmallActivationChange' num2str(i) '_' '12mm''.mat']);
%     for j = 1 : 20
%         
%         
%
%         m = Model{j};
%         c = Config{j};
%         t = SimResults_Time{j};
%         y = SimResults_Y{j};
%         
%         
%        [~,maxAmp(j),~] = c.getOutputOfInterest(t,y);
%         %         m.plot(t,y,'Vid',['Activation' num2str(c.Fr) 'Hz' num2str(c.Amp) 'mm' '_' num2str(j)])
%         %         m.plotActivation;
%         hold on
%         c.getOutputofInterest(t,y);
%         legend('Activation','Mitte','Letzter','Erster','193','355','Location','NorthEast')
%         title(['MaxActivation = ' num2str(1/j^2) 'ActivationTime = ' num2str(5*j)]);
%         hold off
%         print('-f2','-djpeg',['ActivationPlot' num2str(c.Fr) 'Hz' num2str(c.Amp*10)  'mm' '_' num2str(j)]);
%         close all
%         
%     end
%     frequency = linspace(50,15,70);
%     maxAmps = [maxAmps(1:(i-1)*20) maxAmp];
%     maxActivation = linspace(1,.05,20); 
% end