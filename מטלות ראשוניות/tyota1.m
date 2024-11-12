%% introduction to neuronal networks - assignment 1:

Eps = 24;
R1 = 2;
R2 = 0:22;

% in line:

R_t_l = R1 + R2;
I_t_l = Eps./R_t_l;

figure();
hold on;
plot(R2, I_t_l, '.-', MarkerSize = 10, MarkerEdgeColor = 'R');
xlabel('crime scores');
ylabel('number of cities who gave a score in this range');
title('the distribution of crime scores for all US cities');
hold off;



% in parallele:

R_t_p = 1./(1/R1 + 1./R2);
I_t_p = Eps./R_t_p;

figure();
hold on;
plot(R2, I_t_p, '.-', MarkerSize = 10, MarkerEdgeColor = 'R');
xlabel('crime scores');
ylabel('number of cities who gave a score in this range');
title('the distribution of crime scores for all US cities');
hold off;


x = linspace(0,22);
z = 1./(1/2 + 1./x);
y = 24./z;
figure
plot(x, y, 'DisplayName','$C_{\phi}$')
xline(5,'-r', 'DisplayName','$C_{\theta} \rightarrow \infty$')
legend('Location','best', 'Interpreter','latex')