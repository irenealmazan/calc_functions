classdef Phretrieval_functions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [original_object,sup,dp] = prepare_data_ER_HIO(rho,data_exp)
            % this function prepares the data for test 4 conditions: 2DFT dp
            % from original object, no noise and shapr support
            
            original_object = rho;
            sup = zeros(size(rho));%abs(original_object);
           
            range = [-20:20]+65;
            
            range2 = [-20:20]+65;
            
            sup(range,range,range2) = ones(numel(range),numel(range),numel(range2));
            
            
            for ii=1:numel(data_exp)
                dp(:,:,ii) = data_exp(ii).I;
            end
        end

        function [err,support_new] = optimize_support(rho,threshold,FWHM_array,probe,data,angle_list,ki_o,kf_o,X,Y,Z)
            
            for jj = 1:numel(threshold)
                 support_new = Phretrieval_functions.shrink_wrap_support(abs(rho),threshold(jj),FWHM_array,X,Y,Z);
                 err(jj) = DiffractionPatterns.calc_error_multiangle(probe,rho.*support_new,data,angle_list,ki_o,kf_o,X,Y,Z);
            end
            
        end
         
        function [support] = make_support(corners,facet_spacing,edgepad)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            
            support_edge = facet_spacing * edgepad;
            
            %make parallel planes parallel to each pair of facets
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(1,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T1 = (T1>-support_edge/2 & T1<support_edge/2); %two parallel lines
            
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T2 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T2 = (T2>-support_edge/2 & T2<support_edge/2); %two parallel lines
            
            v1 = [corners(4,:) - corners(9,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(9,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T3 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T3 = (T3>-support_edge/2 & T3<support_edge/2); %two parallel lines
            
            v3 = [0 1 0];
            T4 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T4 = (T4>-support_edge/2 & T4<support_edge/2); %two parallel lines
            
            
            support = T1&T2&T3&T4;
            
            %support = T3;
            %support = abs(NW);
            
            
        end
       
         function [new_support] = shrink_wrap_support(mod_object,threshold,FWHM_array,X,Y,Z)
            % this function updates the support by convolving the updated rho
            % with a gaussian function and then selecting the points which have
            % an intensity above a threshold
           
            
           mod_object_blur = Phretrieval_functions.smooth_support(mod_object,FWHM_array,X,Y,Z);
           new_support = (mod_object_blur > threshold*max(mod_object_blur(:)));
            
        end
        
        function [support_new] = smooth_support(support,FWHM_array,X,Y,Z)
            % this function updates the support by convolving the updated rho
            % with a gaussian function and then selecting the points which have
            % an intensity above a threshold
           
            
            [Npix_x,Npix_y,Npix_z] = size(X);
                      
            support_fft = fftshift(fftn(fftshift(support)));
            
            FWHM_x = FWHM_array(1);
            FWHM_y = FWHM_array(2);
            FWHM_z = FWHM_array(3);
            
            gauss_to_convolve = exp(-FWHM_x.*X.^2/Npix_x-FWHM_y.*Y.^2/Npix_y-FWHM_z.*Z.^2/Npix_z);
            
            gauss_fft = fftshift(fftn(fftshift(gauss_to_convolve)));
            
            support_fft_gauss = support_fft.* gauss_fft;
            
            support_new = fftshift(ifftn(fftshift(support_fft_gauss)));
            
            %max_of_int = max(max(max(support_conv_gauss)));
            
            %nonzeros_elements = find(support_conv_gauss> threshold*max_of_int);
            
            %new_support = zeros(Npix_x,Npix_y,Npix_z);
            %new_support(nonzeros_elements) = 1;
            
        end
        
        function [scale_fact] = ini_guess_scalefactor(probe, rho,angle_list_ini, data_exp,ki,kf,X,Y,Z)
            % This function estimates the scaling factor by which the initial guess
            % has to be calculated in order to establish good initial conditions
            % for phase retrieval.
            
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(angle_list_ini,ki,kf,kf-ki);
            
            % modulus square of the wavefront
            [Psi_mod] = DiffractionPatterns.calc_dp( dq_shift,probe,rho,X,Y,Z);
            
            denom = 0;
            numerator = 0;
            for jj = 1:numel(Psi_mod)
               denom = denom + sum(sum(Psi_mod(jj).Psi_mod));
               numerator = numerator + sum(sum(sqrt(Psi_mod(jj).Psi_mod.*data_exp(jj).I)));
            end
            
            scale_fact = numerator/denom;
            
        end
                
        function [rho_new,beta_rho,norm_gradient,gPIEiter,direction_rho] = rho_update(probe, rho,gPIEiter_old,direction_rho_old,angles_list,support,niter, data_exp,depth,err_0,freq_restart,tau_backtrack,beta_ini,counter_max,ki,kf,X,Y,Z)
            % this functions updates rho
            
            % gradient calculation
            [gPIEiter] = GeneralGradient.calc_grad_multiangle(probe, rho,angles_list,data_exp,ki,kf,X,Y,Z);
            
             % calculate the norm of the gradient:
            norm_gradient = gPIEiter(:)'*gPIEiter(:);   
            
            % conjugate parameter calculation
            if niter == 1 || mod(niter,freq_restart) == 0               
                gamma = 0;
            else
               norm_gradient_old = gPIEiter_old(:)'*gPIEiter_old(:);
               gamma = norm_gradient/norm_gradient_old ;
            end
            
            % direction of the scaled steepest descent:
            D = 1/(max(max(max(abs(probe).^2))));
            
            
            %       direction_rho for conjugated gradient ;
            direction_rho = (D/depth)*(-gPIEiter+gamma*direction_rho_old).*support;
            
            % calculate the adaptative step length
            [beta_rho] = GeneralGradient.calc_beta_adaptative_step(probe, rho,angles_list,data_exp,gPIEiter,err_0,direction_rho,tau_backtrack,beta_ini,counter_max,'rho',ki,kf,X,Y,Z);
            
           
            
            % update the object:
             rho_new = rho + beta_rho * direction_rho;
            
             
             
            %{
              % calculate the modulus projection:
            %Pmrho = rho + beta_rho * direction_rho;
             
            % apply constraint of support in real space
%            if ERflag ==1
% 
%                rho_new = support.*Pmrho;
%            else
%                %HIO:
%                rho_new = support.*Pmrho + (1-support).*(rho-0.9*Pmrho);
% 
%            end
             
             %}
            
            
        end
        
        function [dth_new,dq_shift, grad_final_theta,norm_grad_theta,beta] = theta_update(probe, rho,angles_list,data_exp,Niter_theta,error_0,tau_backtrack,beta_ini,counter_max,ki,kf,X,Y,Z,flagDebug)
            %%% this function calculates the gradient of the error metric with
            %%% respect to the position of the angles analytically, and
            %%% the correction theta  step
                        
            qbragg = kf - ki;
                     
            % initialize varibles:
            dth_new = angles_list;
            
            dq_shift = zeros(numel(data_exp),3);
            
            grad_final_theta = zeros(Niter_theta,numel(data_exp));
            
            for ntheta = 1:Niter_theta
                                
                [grad_final_theta(ntheta,:)] = GeneralGradient.calc_grad_theta(probe, rho, data_exp, dth_new,ki,kf,X,Y,Z,flagDebug);
                
                
                for ii = 1:numel(data_exp)%
                    grad_final_theta(ntheta,ii) = grad_final_theta(ntheta,ii)./sum(sum(data_exp(ii).I));
                    if flagDebug == 1
                        display(['dth_new = ' num2str(dth_new(ii)) 'gradient = ' num2str(grad_final_theta(ntheta,ii))]);
                    end
                end
                
                % define the direction of the descent:
                direction = -squeeze(grad_final_theta(ntheta,:)');
                
                % calculate an adaptative step size:
                [beta] = GeneralGradient.calc_beta_adaptative_step(probe, rho, dth_new,data_exp,grad_final_theta,error_0,direction,tau_backtrack,beta_ini,counter_max,'theta',ki,kf,X,Y,Z);

                % corrected theta :
                dth_new = dth_new + beta* direction;
                
                % corrected dqshift:
                [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth_new,ki,kf,qbragg);
                
                % norm of the gradient:
                norm_grad_theta = grad_final_theta(ntheta,:)*grad_final_theta(ntheta,:)';
                
            end
                                    
            
            %%% TEST OF  GRADIENT
            %{
            for jj= dthsearchind(2:end-1)
                grad_manual = test_grad_theta_manually(jj,thscan,fitQerr,data_exp(ii),bmtemp,rho,grad_final_theta(ii).grad(jj));
            end
        
            %}
            
            
        end
        
        function [dth_disp] = generate_angular_jitter(delta_thscanvals,index_to_distort,percentage)
           % this function calculates a distribution of angular gitters
           % between 0 and percentage*diff(delta_thsanvals)
            
           delta_theta = unique(diff(delta_thscanvals));
           delta_theta_ini = delta_theta(1); % corresponds to 100 % jitter
           
           dth_disp = zeros(numel(delta_thscanvals),1);

           
           for jj =1:numel(index_to_distort)
               dth_disp(index_to_distort(jj)) = (delta_theta_ini*percentage/100)*(-1+2*rand);                
           end
            
        end
        
        function  [retrphase,newobj,err_ERHIO] = do_ERHIO(err_ERHIO,dp,support,newobj,er_iter,NW,angles_list,ki,kf,probe,d2_bragg,X,Y,Z,plotResults,flagER,flagER_direct)
            
            
            
            for jj = 1:er_iter
                if flagER
                    [retrphase newobj] = erred3( sqrt(dp),support,1, 10000, newobj,plotResults);
                    string_iter = ['ER_iter'];
                else
                    [retrphase newobj] = hio3( sqrt(dp),support,1, 10000, newobj,0.7,plotResults);
                    string_iter = ['HIO_iter'];
                end
                
                if flagER_direct
                    
                    finalobj_3DFT = (ifftn(newobj.dp));
                    
                    finalobj_2DFT = DiffractionPatterns.From3DFT_to_2DFT(finalobj_3DFT,angles_list,probe,ki,kf,X,Y,Z);
                    support_2DFT = DiffractionPatterns.From3DFT_to_2DFT(support,angles_list,probe,ki,kf,X,Y,Z);
                    
                    [err,index_NW,val_NW,index_rho,val_rho] = Phretrieval_functions.decide_flip(NW,finalobj_2DFT,support_2DFT,angles_list,ki,kf,d2_bragg,X,Y,Z);
                  
                    
                else
                    err = 0;
                end
                
                err_ERHIO = [err_ERHIO err];
                
                display([string_iter num2str(numel(newobj.chi)) ' error: ' num2str(err) ' chi value: ' num2str(newobj.chi(end)) ' \n'])
                
               
            end
            
            
        end
               
        function [err,index_NW,val_NW,index_rho,val_rho] = decide_flip(NW,rho,support,angles_list,ki,kf,d2_bragg,X,Y,Z)
            % this function checks wether the retrieved object is flipped
            % with respect to the original one.
            
            % conjugate the diffraction pattern to flip the object and then
            % shift
            finalobj = (ifftn(conj(fftn(rho))));
            [finalobj_shift] = DiffractionPatterns.shift_object(NW,finalobj,angles_list,ki,kf,kf-ki,d2_bragg,X,Y,Z);       
            err_1 = DiffractionPatterns.calculate_error_realspace(abs(NW),abs(finalobj_shift));

            % only shift the object
            finalobj_2 = ifftn(fftn(rho));
            [finalobj_2_shift] = DiffractionPatterns.shift_object(NW,finalobj_2,angles_list,ki,kf,kf-ki,d2_bragg,X,Y,Z);            
            err_2 = DiffractionPatterns.calculate_error_realspace(abs(NW),abs(finalobj_2_shift));
            
            if err_1<err_2
                finalobj =finalobj_shift;
                support_final = (ifftn(conj(fftn(support))));
            else
                finalobj =  finalobj_2_shift;
                support_final = support;
            end
            
            support_shift = DiffractionPatterns.shift_object(NW,support_final,angles_list,ki,kf,kf-ki,d2_bragg,X,Y,Z); 
            support_shift_abs = abs(support_shift);            
            support_shift_final = abs(support_shift>0.1*max(abs(support_shift(:)))) ;
            % phase of the final object at the center pixel:
            Nx_c = round(size(finalobj,1)/2);
            Ny_c = round(size(finalobj,2)/2);
            Nz_c = round(size(finalobj,1)/2);
            phase_offset_finalobj = angle(finalobj(Nx_c,Ny_c,Nz_c));
            
             % phase of the real object at the center pixel:
            phase_offset_NW = angle(NW(Nx_c,Ny_c,Nz_c));
            
            
            %err = DiffractionPatterns.calculate_error_realspace(NW*exp(-1i*phase_offset_NW),finalobj*exp(-1i*phase_offset_finalobj).*support_shift_final);
            err = DiffractionPatterns.calculate_error_realspace(NW*exp(-1i*phase_offset_NW),finalobj*exp(-1i*phase_offset_finalobj));
            
            % save non-zero values           
            [index_rho] = find(abs(finalobj(:).*support_shift(:))>1e-16);
            val_rho = finalobj(index_rho);
            index_NW = index_rho;
            val_NW = NW(index_NW);
            
        end
        
        
    end
        
   
        
        
end
    

