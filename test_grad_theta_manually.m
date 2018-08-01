function grad_manual = test_grad_theta_manually(mmtest,dth_grid,err_single,data_exp,beam,rho,grad_theta_calc)

    % This function computes the gradient of the error with respect to theta manually and using the function
    % calc_grad_theta
    % and plots the result

    global  ki_o kf_o X Y Z

    % theta + delta_theta
    
     dth =  dth_grid(mmtest+1) ;
%     
     Ry = [cosd(-dth) 0 sind(-dth);
        0 1 0;
        -sind(-dth) 0 cosd(-dth)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';

    qbragg = kf_o - ki_o;
    
    dqtest_up = kf - ki - qbragg;


     Qterm_up = exp(1i* dqtest_up(1) * X) .* ...
            exp(1i* dqtest_up(2) * Y) .* ...
            exp(1i* dqtest_up(3) * Z);
    
    
    [err_single_plus] = calc_error_theta_singlepos(rho,beam,data_exp,Qterm_up);
    
    % theta + delta-theta
    
     dth =  dth_grid(mmtest-1) ;
%     
     Ry = [cosd(-dth) 0 sind(-dth);
        0 1 0;
        -sind(-dth) 0 cosd(-dth)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';

    qbragg = kf_o - ki_o;
    
    dqtest_minus = kf - ki - qbragg;


     Qterm_minus = exp(1i* dqtest_minus(1) * X) .* ...
            exp(1i* dqtest_minus(2) * Y) .* ...
            exp(1i* dqtest_minus(3) * Z);
    
    [err_single_minus] = calc_error_theta_singlepos(rho,beam,data_exp,Qterm_minus);
    
    % rough derivative:
    
    grad_manual = (err_single_plus - err_single_minus)/(dth_grid(mmtest+1)-dth_grid(mmtest-1));
    
    
    
    
    % plot derivative:
   
    figure(27);clf;
    hold on;
    subplot(2,1,1);
    plot(dth_grid,err_single,'LineWidth',3.0);
    hold on;
    cst1 = err_single(mmtest)-grad_manual*dth_grid(mmtest);
    plot(dth_grid,grad_manual.*dth_grid+cst1);
    plot(dth_grid(mmtest+1),err_single_plus,'or');
    plot(dth_grid(mmtest-1),err_single_minus,'or');
    
    subplot(2,1,2);
    plot(dth_grid,err_single,'LineWidth',3.0);
    hold on;
    cst2 = err_single(mmtest)-grad_theta_calc*dth_grid(mmtest);
    plot(dth_grid,+grad_theta_calc.*dth_grid+cst2);

end