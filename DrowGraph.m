function DrowGraph(result)
    set(groot, 'DefaultTextInterpreter', 'Latex');
    set(groot, 'DefaultLegendInterpreter', 'Latex');
    figure(1);
    title('Self-Localization', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value', 'UKF','QS-ARUKF', 'The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    grid on;
    axis equal;
    
    figure(2);
    title('UKF', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value','UKF', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    grid on;
    axis equal;
    
    figure(3);
    title('QS-ARUKF', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value','QS-ARUKF', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    grid on;
    axis equal;
    
    figure(4);
    title('The Proposed method', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value','The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    grid on;
    axis equal;
    
    figure(5);
    title('Self-Localization', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value','GPS','UKF','QS-ARUKF', 'The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    xlim([-20 70]); 
    ylim([-8 70]);
    grid on;
    axis equal;
    
    figure(6);
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(result.xTrue(:, 1), result.xTrue(:, 2),'k','linewidth', 1.5); hold on;
    plot(result.z(:, 1), result.z(:, 2), '.r'); hold on;
    plot(result.xEst_UKF(:, 1), result.xEst_UKF(:, 2),'.b','linewidth', 1.5); hold on;
    plot(result.xEst_AMUKF(:, 1), result.xEst_AMUKF(:, 2),'.g','linewidth', 1.5); hold on;
    plot(result.xEst_QS_ARUKF(:, 1), result.xEst_QS_ARUKF(:, 2),'.m','linewidth', 1.5); hold on;
    title('Self-Localization', 'fontsize', 16, 'interpreter','latex');
    xlabel('X [m]','fontsize', 16,'interpreter','latex')
    ylabel('Y [m]','fontsize', 16,'interpreter','latex')
    leg=legend('True value','GPS','UKF','QS-ARUKF','The Proposed method','Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    grid on;
    axis equal;
    
    figure(7);
    subplot(3,1,1);
    plot(result.time,result.xTrue(:, 1),'k',result.time,result.xEst_UKF(:, 1),'b--', ...
         result.time,result.xEst_AMUKF(:, 1),'g--', result.time,result.xEst_QS_ARUKF(:, 1),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('X[m]','fontsize', 12,'interpreter','latex')
    xlim([0 60]);
    grid on
    subplot(3,1,2);
    plot(result.time,result.xTrue(:, 2),'k',result.time,result.xEst_UKF(:, 2),'b--', ...
         result.time,result.xEst_AMUKF(:, 2),'g--', result.time,result.xEst_QS_ARUKF(:, 2),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('Y[m]','fontsize', 12,'interpreter','latex')
    xlim([0 60]);
    grid on
    subplot(3,1,3);
    plot(result.time,result.xTrue(:, 3),'k',result.time,result.xEst_UKF(:, 3),'b--', ...
         result.time,result.xEst_AMUKF(:, 3),'g--', result.time,result.xEst_QS_ARUKF(:, 3),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    xlabel('Time [s]','fontsize', 12,'interpreter','latex')
    ylabel('$\theta$[rad]','fontsize', 12,'interpreter','latex')
    leg=legend('True value', 'UKF','QS-ARUKF', 'The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    xlim([0 60]);
    grid on

    figure(8);
    subplot(3,1,1);
    plot(result.time,result.RMSE_UKF(:, 1),'b--',result.time,result.RMSE_AMUKF(:, 1),'g--', result.time,result.RMSE_QS_ARUKF(:, 1),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('Position Error X[m]','fontsize', 10,'interpreter','latex')
    leg=legend('UKF','QS-ARUKF', 'The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    xlim([0 60]);
    grid on
    subplot(3,1,2);
    plot(result.time,result.RMSE_UKF(:, 2),'b--',result.time,result.RMSE_AMUKF(:, 2),'g--', result.time,result.RMSE_QS_ARUKF(:, 2),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('Position Error Y[m]','fontsize', 10,'interpreter','latex')
    xlim([0 60]);
    grid on
    subplot(3,1,3);
    plot(result.time,result.RMSE_UKF(:, 3),'b--',result.time,result.RMSE_AMUKF(:, 3),'g--',result.time,result.RMSE_QS_ARUKF(:, 3),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    xlabel('Time [s]','fontsize', 12,'interpreter','latex')
    ylabel('Posture Error $\theta$[rad]','fontsize', 10,'interpreter','latex')
    xlim([0 60]);
    grid on
    
    figure(9);
    subplot(3,1,1);
    plot(result.time,result.RMSE_UKF(:, 1),'b--',result.time,result.RMSE_AMUKF(:, 1),'g--', result.time,result.RMSE_QS_ARUKF(:, 1),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('Position Error X[m]','fontsize', 10,'interpreter','latex')
    leg=legend('UKF','ARUKF', 'The Proposed method', 'Location','best','fontsize', 12);
    set(leg,'interpreter','latex');
    xlim([28 38]);
    ylim([-5 5]);
    grid on
    subplot(3,1,2);
    plot(result.time,result.RMSE_UKF(:, 2),'b--',result.time,result.RMSE_AMUKF(:, 2),'g--', result.time,result.RMSE_QS_ARUKF(:, 2),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    ylabel('Position Error Y[m]','fontsize', 10,'interpreter','latex')
    xlim([28 38]);
    ylim([-5 5]);
    grid on
    subplot(3,1,3);
    plot(result.time,result.RMSE_UKF(:, 3),'b--',result.time,result.RMSE_AMUKF(:, 3),'g--',result.time,result.RMSE_QS_ARUKF(:, 3),'m--');
    set(gca, 'fontsize', 10, 'fontname', 'times');
    xlabel('Time [s]','fontsize', 12,'interpreter','latex')
    ylabel('Posture Error $\theta$[rad]','fontsize', 10,'interpreter','latex')
    xlim([28 38]);
    ylim([-5 5]);
    grid on
    
    RMSE = @(x) sqrt(mean(x.^2));
    fprintf('%10s %10s\n','variable','RMSE(ekf)');
    for p=1:3
        vname = sprintf('x%d',p);
        fprintf('%10s %10.5f \n',vname,RMSE(result.xTrue(:, p) - result.xEst_UKF(:, p)));
        fprintf('%10s %10.5f \n',vname,RMSE(result.xTrue(:, p) - result.xEst_AMUKF(:, p)));
        fprintf('%10s %10.5f \n',vname,RMSE(result.xTrue(:, p) - result.xEst_QS_ARUKF(:, p)));
    end
end