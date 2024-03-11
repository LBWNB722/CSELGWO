clc;
clear all;
n=50;
Function_name='F1'; 
maxiter=2000;
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Best_score,Best_pos,Convergence_curve]= CSELGWO(n,maxiter,lb,ub,dim,fobj);
figure('Position',[500 500 660 290])
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('CSELGWO')

display(['The best solution obtained by CSELGWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by CSELGWO is : ', num2str(Best_score)]);