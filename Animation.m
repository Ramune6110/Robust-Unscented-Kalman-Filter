function Animation(UKF, AMUKF, QS_ARUKF, i)
    figure(1);
    if i >= 250 && i <= 400
        plot(UKF.xTrue(1), UKF.xTrue(2), '.k','linewidth', 1.5); hold on;
        plot(UKF.xEst(1), UKF.xEst(2), '.b','linewidth', 1.5); hold on;
        plot(AMUKF.xEst(1), AMUKF.xEst(2), '.g','linewidth', 1.5); hold on;
        plot(QS_ARUKF.xEst(1), QS_ARUKF.xEst(2), '.m','linewidth', 1.5); hold on;
        if rem(i, 20) == 0
            ShowErrorEllipse(UKF, 0);
            axis equal;
            grid on;
            drawnow;
        end
        if rem(i, 40) == 0
            ShowErrorEllipse(AMUKF, 1);
            ShowErrorEllipse(QS_ARUKF, 2);
            axis equal;
            grid on;
            drawnow;
        end
    end
    
    figure(2);
    if i >= 200 && i <= 450
        plot(UKF.xTrue(1), UKF.xTrue(2), '.k','linewidth', 1.5); hold on;
        plot(UKF.xEst(1), UKF.xEst(2), '.b','linewidth', 1.5); hold on;
        if rem(i, 5) == 0
            ShowErrorEllipse(UKF, 0);
            axis equal;
            grid on;
            drawnow;
        end
    end
    
    figure(3);
    if i >= 200 && i <= 450
        plot(UKF.xTrue(1), UKF.xTrue(2), '.k','linewidth', 1.5); hold on;
        plot(AMUKF.xEst(1), AMUKF.xEst(2), '.g','linewidth', 1.5); hold on;
        if rem(i, 30) == 0
            ShowErrorEllipse(AMUKF, 1);
            axis equal;
            grid on;
            drawnow;
        end
    end
    
    figure(4);
    if i >= 200 && i <= 450
        plot(UKF.xTrue(1), UKF.xTrue(2), '.k','linewidth', 1.5); hold on;
        plot(QS_ARUKF.xEst(1), QS_ARUKF.xEst(2), '.m','linewidth', 1.5); hold on;
        if rem(i, 30) == 0
            ShowErrorEllipse(QS_ARUKF, 2);
            axis equal;
            grid on;
            drawnow;
        end
    end
    
    figure(5);
    if i >= 200 && i <= 400
        plot(UKF.xTrue(1), UKF.xTrue(2), '.k','linewidth', 1.5); hold on;
        plot(UKF.xEst(1), UKF.xEst(2), '.b','linewidth', 1.5); hold on;
        if i >= 300 && i <= 310
        plot(UKF.z(1), UKF.z(2), '.r','linewidth', 1.5); hold on;
        end
        plot(AMUKF.xEst(1), AMUKF.xEst(2), '.g','linewidth', 1.5); hold on;
        plot(QS_ARUKF.xEst(1), QS_ARUKF.xEst(2), '.m','linewidth', 1.5); hold on;
        axis equal;
        grid on;
        drawnow;
        if i >= 290 && i <= 330
            if rem(i, 10) == 0
                ShowErrorEllipse(UKF, 0);
                axis equal;
                grid on;
                drawnow;
            end
            if rem(i, 20) == 0
                ShowErrorEllipse(AMUKF, 1);
                ShowErrorEllipse(QS_ARUKF, 2);
                axis equal;
                grid on;
                drawnow;
            end
        end
    end
end