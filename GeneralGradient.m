classdef GeneralGradient
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
             
        function [gradtot] = calc_grad_multiangle(probe, rho,angles_list, data,ki_ini,kf_ini,X,Y,Z)
            
            gradtot = zeros(size(rho));
            
            [dqshift] = DiffractionPatterns.calc_dqshift_for_given_th(angles_list,ki_ini,kf_ini,kf_ini-ki_ini);
            
            [Psi_mod,~,~,Psij,Qterm] = DiffractionPatterns.calc_dp(dqshift,probe,rho,X,Y,Z);
            
            for ii=1:numel(data)
 
                Psig = sqrt(data(ii).I).*exp(1i*angle(Psij(ii).FT_Psij));%
                grad = Psij(ii).FT_Psij - Psig; %%%% NEW % Psij - Psig;
                grad = fftshift(ifftn(fftshift(grad)));
                grad = repmat(grad, [1 1 size(probe,3)]);
                grad = grad .* conj(Qterm(ii).Qterm);
                grad = grad .* conj(probe);
                grad = 2*grad; %add factor of 2 to make sure that gradient is equivalent
                
                gradtot = gradtot + grad;
            end
        end
    
        function [grad_final_theta] = calc_grad_theta(probe, rho, data, dth_nominal,ki,kf,X,Y,Z,flagDebug)
            % this function calculates the gradient of the error metric with
            % respect to theta. dth_delta = 0 if we are using the function to
            % refine the angular position.    
            
            qbragg = kf - ki;
            
            % calculate the current momentum transfer:
            [dq_shift,~,~,dq_shift_deriv] = DiffractionPatterns.calc_dqshift_for_given_th(dth_nominal,ki,kf,qbragg);
                     
            % calculate the derivative of dq_shift with respect to theta:
            %[dq_shift_deriv] = DiffractionPatterns.calc_dq_deriv_analytical(dth_nominal,ki,kf,qbragg);
             
            % display the derivative of dq_shift
            if flagDebug == 1
                for jj = [1:round(numel(dth_nominal)/5):numel(dth_nominal)]
                    h2(jj).h = DisplayFunctions.display_calc_dqshift_deriv(dth_nominal(jj),dq_shift_deriv(jj,:),dq_shift(jj,:),ki,kf,30+jj,'analitically');
                end
            end
            [Psi_mod,~,~,FT_Psij,Qterm] = DiffractionPatterns.calc_dp(dq_shift,probe,rho,X,Y,Z);
            
            grad_final_theta = zeros(numel(dth_nominal),1);
           
            for jj = 1:numel(dth_nominal)
                
                % derivative of Qterm with respect to dq_shift
                deriv_Qterm_theta = (dq_shift_deriv(jj,1).*X + dq_shift_deriv(jj,3).*Z).*Qterm(jj).Qterm;
                
                % derivative of the diffraction pattern with respect to dq_shift
                Psij_deriv_theta = probe.*rho.*deriv_Qterm_theta;
                
                Psij_deriv_theta = sum(Psij_deriv_theta,3);
                
                Psij_deriv_theta = fftshift(fftn(fftshift(Psij_deriv_theta)));
                
                Im_Psij_deriv_theta = imag(conj(FT_Psij(jj).FT_Psij).*Psij_deriv_theta);
                
                % calculate the gradient:
                
                grad_theta = (1 -sqrt(data(jj).I./(Psi_mod(jj).I))).* Im_Psij_deriv_theta; %%%% NEW
                
                grad_final_theta(jj) = -2*sum(sum(grad_theta))/(numel(grad_theta));
                
            end
            
            %%% check that you calculate the gradient:           
            if 0
                for jj = [1:round(numel(dth_nominal)/5):numel(dth_nominal)]
                    error_single(jj) = DiffractionPatterns.calc_error_multiangle(probe, rho, data(jj), dth_nominal(jj),ki,kf,X,Y,Z);
                    DisplayFunctions.display_calc_grad_theta(probe, rho, data(jj), dth_nominal(jj),grad_final_theta(jj),error_single(jj) ,X,Y,Z,ki,kf,40+jj);
                end
            end
            %        
            %display_calc_grad2_theta(probe, rho, data, dth_nominal,grad2_final_theta,grad_final_theta)
        end
           
        function [beta_iter,err_plusalpha,beta_track] = calc_beta_adaptative_step(probe, rho,angles_list,data,gradtot,err_0,direction,tau_backtrack,beta_ini,counter_max,flag,ki,kf,X,Y,Z)
            % this function calculates the adaptative step length for the
            % correction of the rho. We folllow Nocedal "Backtracking Line Search
            % The imputs are the following:
            %   probe: information about the beam (important in ptychography)
            %   rho: the object in its current state (before the update)
            %   data: structure containing the angles
            %   gradtot: vector with all the gradients for each angle
            
            
            
            % Armijo's condition: if err(alpha_i) > err(alpha=0) +
            % c1*deriv_err_with_theta *alpha_i
            
            % calculate the value of alpha_ini that we choose as the value for
            % which the linear approximation of the error metric becomes positive
            c1 = 1e-3; % see Nocedal
            counter = 1; % counter to track the evolution of the error, alpha and the approximation of the error with the iterations            
            err_linear_aprox(counter) = - err_0;
            
            slope_alpha_0 = c1*real(direction(:)'*gradtot(:));
            
            
            % rho/theta at initial beta           
            switch flag
                
                case 'rho'
                    rho_beta = rho + beta_ini * direction;
                    theta_beta = angles_list;                    
                case 'theta'                    
                    rho_beta = rho;
                    theta_beta = angles_list + beta_ini * direction;
            end
            
            % estimate the value of the error at rho + alpha_ini*gradtot_rho
            err_plusalpha(counter) =  DiffractionPatterns.calc_error_multiangle(probe, rho_beta, data,theta_beta,ki,kf,X,Y,Z);
            
            Delta_err(counter) = err_plusalpha(counter) - err_0;
            
            beta_iter = beta_ini;
            
            beta_track(counter) = beta_iter;
            
            while( Delta_err(counter) >= beta_iter*slope_alpha_0)
                
                % update the counter
                counter = counter + 1;
                
                % update beta
                beta_iter = beta_iter*tau_backtrack;
                
                % move rho or theta
                switch flag                 
                    case 'rho'
                        rho_beta = rho + beta_iter * direction;                
                    case 'theta'
                        theta_beta = angles_list + beta_iter*direction;
                end
                
                
                % estimate the value of the error at rho + alpha_ini*gradtot_rho
                err_plusalpha(counter) =  DiffractionPatterns.calc_error_multiangle(probe, rho_beta, data,theta_beta,ki,kf,X,Y,Z);
                
                
                %recalculate the difference between error metric at alpha_iter and
                %error metric at alpha = 0
                Delta_err(counter) = err_plusalpha(counter) - err_0;
                
                % calculate the linear approximation of the error metric
                err_linear_aprox(counter) =  beta_iter*slope_alpha_0;
                
                display(['err(' num2str(beta_iter) ' ) - err(0) = ' num2str(Delta_err(counter)) ' and linear approx. ' num2str(err_linear_aprox(counter))])
                
                beta_track(counter) = beta_iter;
                
                if counter > counter_max
                    display('beta not found');
                    break;
                    
                end
                
            end
            
            
            
            
        end
        
        function [err_plusalpha] = calc_err_vs_rho_theta(probe, rho,angles_list,data,direction_rho,direction_theta,ki,kf,X,Y,Z)
            % this function calculates the contour plot of the error metric
            % as a function of rho and theta
            
            beta_rho = [1e-3 1e-2 1e-1 1];
            beta_theta = [1e-3 1e-2 1e-1 1]*1e-5;
            
            err_plusalpha = zeros(numel(beta_rho),numel(beta_theta));
            
            for kk = 1:numel(beta_rho)
                rho_beta = rho + beta_rho(kk) * direction_rho;                
                for jj = 1:numel(beta_theta)                    
                    theta_beta = angles_list + beta_theta(jj)*direction_theta;
                    % estimate the value of the error at rho + alpha_ini*gradtot_rho
                    err_plusalpha(kk,jj) =  DiffractionPatterns.calc_error_multiangle(probe, rho_beta, data,theta_beta,ki,kf,X,Y,Z);                    
                end                
            end            
        end
        
        
    end
        
   
        
        
end