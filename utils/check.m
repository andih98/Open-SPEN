function [rep] = check(seq,TE,TR,delayTR)
    [ok, error_report]=seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    try 
    % exportgraphics((seq.plot('TimeRange', [0 (TR-delayTR)])),'seq.png','Resolution',1000)
    seq.plot('stacked',1,'timeRange', [0 TR-delayTR]);
    catch
    % exportgraphics((seq.plot('TimeRange', [0 (TR-delayTR)])),'seq.png','Resolution',1000)
    seq.plot('stacked', 1);
    end
    
    % calculate trajectory 
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
    
    % plot k-spaces
    figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
    hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
    
    % Vorhandener X-Y Plot
    figure; plot(ktraj(1,:),ktraj(2,:),'b',...
                 ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
    title('k-Space Trajectory (X-Y Plane)');
    xlabel('k_x [1/m]');
    ylabel('k_y [1/m]');
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    grid on;

    % --- HIER ERWEITERT: X-Z Plot ---
    % Plot der X-Z-Ebene (kx vs kz)
    figure;
    hold on;
    % Gesamte Trajektorie (Linie)
    plot(ktraj(1, :), ktraj(3, :), 'b');
    % ADC-Abtastpunkte (Punkte)
    plot(ktraj_adc(1, :), ktraj_adc(3, :), 'r.', 'MarkerSize', 6);
    hold off;
    
    title('k-Space Trajectory (X-Z Plane)');
    xlabel('k_x [1/m]');
    ylabel('k_z [1/m]');
    legend('Trajectory', 'ADC-Samples', 'Location', 'best');
    axis('equal'); % Wichtig für korrekte Proportionen
    grid on;
    % --- Ende der Erweiterung ---

    % figure, plot(t_ktraj,ktraj(1,:));
    % sanity checks
    TE_check=(t_refocusing(1)-t_excitation(1))*2;
    fprintf('intended TE=%.03f ms, actual spin echo TE=%.03fms\n', TE*1e3, TE_check*1e3); 
    rep = seq.testReport;
end