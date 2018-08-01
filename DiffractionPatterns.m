classdef DiffractionPatterns
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
       function [dqshift,ki,kf,dq_shift_deriv] = calc_dqshift_for_given_th(dth,ki_ini,kf_ini,qbragg)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            ki = zeros(numel(dth),3); 
            kf = zeros(numel(dth),3); 
            dqshift = zeros(numel(dth),3);
            
            for ii = 1:numel(dth)
                
                [Ry,Ry_deriv] = RotationMatrix.rock_curve(dth(ii),0.0);
                
                ki(ii,:) = (Ry * ki_ini.').';
                kf(ii,:) = (Ry * kf_ini.').';
                
                dqshift(ii,:) = (kf(ii,:)-ki(ii,:))-qbragg;
                dq_shift_deriv(ii,:) = (Ry_deriv * qbragg');

            end
            
        end
        
       function [dq_shift_deriv] = calc_dq_deriv_manually(dth,ki_ini,kf_ini,qbragg)
            % This function calculates the first derivative of the vector dq
            % (linking the different positions of the detector at different theta
            % angles - differents points in the rocking curve)
            
            dq_shift_deriv = zeros(numel(dth),3);
            
            dth_infinitesimal = 1e-4;
            
            for jj = 1:numel(dth)
                dq_shift_plus = DiffractionPatterns.calc_dqshift_for_given_th(dth(jj)+ dth_infinitesimal,ki_ini,kf_ini,qbragg);
                dq_shift_minus = DiffractionPatterns.calc_dqshift_for_given_th(dth(jj)- dth_infinitesimal,ki_ini,kf_ini,qbragg);

                dq_shift_deriv(jj,:) = (dq_shift_plus - dq_shift_minus)/(2*dth_infinitesimal);
            end
            
        end
         
       function [Psi_mod,rock_curve,Psij,FT_Psij,Qterm] = calc_dp(dq_shift,probe,rho,X,Y,Z)
            % This function calculates a diffraction pattern corresponding
            % to one single point of the rocking curve and one single dq
            % shift
           
            rock_curve = zeros(size(dq_shift,1),1);
            
            for jj = 1:size(dq_shift,1)
                % calculate the Q term which depends on the theta angle distorted!!! :
                Qterm(jj).Qterm = exp(1i* dq_shift(jj,1) * X) .* ...
                    exp(1i* dq_shift(jj,2) * Y) .* ...
                    exp(1i* dq_shift(jj,3) * Z);
                
                %%{
                % R(Q*P_j*rho): projected volume
                Psij(jj).Psij = sum( rho.*probe.*Qterm(jj).Qterm, 3);
                
                % 2D_FT: diffracted wave
                FT_Psij(jj).FT_Psij = fftshift(fftn(fftshift( Psij(jj).Psij)));
                
                % modulus square of the diffracted wave
                Psi_mod(jj).I =  FT_Psij(jj).FT_Psij.* conj(FT_Psij(jj).FT_Psij);
                
                % integrated intensity:
                rock_curve(jj) = sum(sum(Psi_mod(jj).I));
                    %}
            end
            
        end
                
       function [errtot] = calc_error_multiangle(probe, rho, data,angle_list,ki,kf,X,Y,Z)
            % This function calculates the diffracted intensities
            % differences for the set of data in the the structure "data"
            
            qbragg = kf - ki;
            
            errtot=0;
            
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(angle_list,ki,kf,qbragg);
            
            [Psi_mod] = DiffractionPatterns.calc_dp( dq_shift,probe,rho,X,Y,Z);
            
            for ii=1:numel(data)
                err = sqrt(data(ii).I) - sqrt(Psi_mod(ii).I);
                err = sum(err(:).^2)/numel(err);
                errtot = errtot + err;
            end
            
            errtot = errtot/numel(data);
            
       end
        
       function [err] = calc_error_3DFT(dp,rho)
           
           A = fftshift(fftn(fftshift(rho)));
           
           %we apply the phase guess to the measured dp amplitudes
           B = abs(dp) .* exp(sqrt(-1) * angle(A));
           
           %length(find(angle(A) ~= angle(B)))
           
           %we bring the new guess to real space via IFT
           C = ifftn(B);
                      
                    
           %return to dp space and apply the
           A = fftn(C);
           
          err = sum(sum(sum( (abs(A)-abs(dp)).^2))) / numel(dp) ;

           
       end
       
        function [Pmrho,Psi_mod] = InvFT_2D(rho,dth,probe,ki_o,kf_o,X,Y,Z)
            % This function makes the inverse Fourier transform (from the
            % reciprocal space to the real space, for each slice of the imput
            
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth,ki_o,kf_o,kf_o-ki_o);
            [Psi_mod,~,~,FT_Psij,Qterm] = DiffractionPatterns.calc_dp(dq_shift,probe,rho,X,Y,Z);
            
            Pmrho = zeros(size(rho));
            
            for jj=1:numel(Psi_mod)
                                
                Pmrho_dummy = fftshift(ifftn(fftshift(FT_Psij(jj).FT_Psij)));
                Pmrho_dummy2 = repmat(Pmrho_dummy,[1 1 numel(Psi_mod)])/numel(Psi_mod);
                Pmrho =  Pmrho + conj(Qterm(jj).Qterm).*Pmrho_dummy2;
                                
            end
        end
        
        function [rho_new] = From3DFT_to_2DFT(rho,dth,probe,ki_o,kf_o,X,Y,Z)
            % This function takes the object calculated in the squeeued
            % frame of the 3DFT to the orthogonal frame of the 2DFT
            FT_rho = fftshift((fftn(fftshift(rho)))); 
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth,ki_o,kf_o,kf_o-ki_o);
            [~,~,~,~,Qterm] = DiffractionPatterns.calc_dp(dq_shift,probe,rho,X,Y,Z);
            
            rho_new = zeros(size(rho));
            
            for jj=1:size(FT_rho,3)
                FT_Psij(jj).FT_Psij = FT_rho(:,:,jj);              
                Pmrho_dummy = fftshift(ifftn(fftshift(FT_Psij(jj).FT_Psij)));
                Pmrho_dummy2 = repmat(Pmrho_dummy,[1 1 size(FT_rho,3)])/size(FT_rho,3);
                rho_new =  rho_new + conj(Qterm(jj).Qterm).*Pmrho_dummy2;
                                
            end
        end
        
        function [FT_rho,dp] = From2DFT_to_3DFT(rho,data)
            % This function takes the object calculated in the squeeued
            % frame of the 3DFT to the orthogonal frame of the 2DFT
            FT_rho = ((fftn((rho))));
            
            for jj = 1:numel(data)
               dp(:,:,jj) = data(jj).I; 
            end
            
           
        end
        
        
        function [Psi_mod,rock_curve] = calc_rock_curve_2DFT(rho,probe,dth,ki_ini,kf_ini,qbragg,X,Y,Z)
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth,ki_ini,kf_ini,qbragg);
            
            [Psi_mod,rock_curve,~,~,~] = DiffractionPatterns.calc_dp(dq_shift,probe,rho,X,Y,Z);
            
        end
        
        function [rock_curve_3D,intens] = calc_rock_curve_3DFT(img,addNWstrain,mncntrate)
            
            if addNWstrain
                fk = fftshift( fftn( fftshift( img ) ));
            else
                fk = fftshift( fftn( fftshift( abs(img) ) ));
            end
            
            intens = fk.*conj(fk);
                        
            midslice = round(size(intens,3)/2);
            
            intens_mean = squeeze(intens(:,:,midslice));
            
            mn_3D = mean(intens_mean(:));
            
            rock_curve_3D = zeros(size(intens,3),1);
            for jj = 1:size(intens,3)
                %rock_curve_3D(jj) = sum(sum(squeeze(intens(:,:,jj)./mn_3D * mncntrate))) ;
                rock_curve_3D(jj) = sum(sum(squeeze(intens(:,:,jj)))) ;
            end
            
        end
      
        function [img_comp] = calc_compatible_rho(img,addNWstrain)
            
            if addNWstrain
                fk = fftshift( fftn( fftshift( img ) ));
            else
                fk = fftshift( fftn( fftshift( abs(img) ) ));
            end
            
            intens = fk.*conj(fk);
            %patt = fftshift( fftn( fftshift( intens ) ) );
            
            mask = find(intens<1);
            fk_mask = fk;
            fk_mask(mask) = 0;
            img_comp = fftshift(ifftn(fftshift(fk_mask)));
            
        end
      
        function [X_recip,Y_recip,Z_recip, qz_pixel_size] = calc_reciprocal_space(dqshift,Npix,depth,d2_bragg)
                       
            % pixel size on the z direction of the reciprocal space:
            qz_pixel_size = abs((dqshift(1,3) - dqshift(end,3))/size(dqshift,1));
            
%             
%              [X_recip,Y_recip,Z_recip] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
%                 [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg),...
%                 [-depth/2:depth/2-1].*2*qz_pixel_size); 
%             
             [X_recip,Y_recip,Z_recip] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg),...
                [-depth/2:depth/2-1].*qz_pixel_size); 
        end
        
        function [autocorr,rho_final_shift,shift_directspace] = calc_autocorr_shift_object(rho_1,rho_2,X,Y,Z,X_recip,Y_recip,Z_recip)
            
            fft_rho_original = fftn(rho_1);
            fft_rho_final = fftn(rho_2);
            
            fft_autocorr = fft_rho_original.*conj(fft_rho_final);
            autocorr = fftshift(ifftn(fft_autocorr));
            
            
              
             % maximum value of the autocorrelation in direct space
            [mxv,idx] = max(autocorr(:));
            [r,c,p] = ind2sub(size(autocorr),idx);
            
            % multiply by phase ramp FT of object:
            shift_directspace = [X(r,c,p) - X(65,65,65) Y(r,c,p) - Y(65,65,65) Z(r,c,p) - Z(65,65,65)];            
            %rho_final_shift = rho_2 + shift_directspace(1);
            fft_rho_final_shift = fft_rho_final.*exp(-1i*(shift_directspace(1)*X_recip + shift_directspace(2)*Y_recip + shift_directspace(3)*Z_recip));

            % use pixels:
            shift_pixels = [r-65 c-65 p-65];
            rho_final_shift_pixel = circshift(rho_2,shift_pixels); 
            
            rho_final_shift = ifftn(fft_rho_final_shift);
            
        end
        
        function [rho_final_shift,shift_directspace] = shift_object(rho_original,rho_final,angles_list,ki,kf,qbragg,d2_bragg,X,Y,Z)
           
             Npix = size(X,1);
            depth = size(X,3);
            
             % calculate a list of momentum transfers:
            [dqshift] = DiffractionPatterns.calc_dqshift_for_given_th(angles_list,ki,kf,qbragg);
            
            [X_recip,Y_recip,Z_recip, qz_pixel_size] = DiffractionPatterns.calc_reciprocal_space(dqshift,Npix,depth,d2_bragg);
           
            [autocorr,rho_final_shift,shift_directspace] = DiffractionPatterns.calc_autocorr_shift_object(rho_original,rho_final,X,Y,Z,X_recip,Y_recip,Z_recip);
           
            
        end
        
        function [object_final_shift] = shift_object_known_shift(object,shift,angles_list,ki,kf,qbragg,X,d2_bragg)
            
            Npix = size(X,1);
            depth = size(X,3);
            
            [dqshift] = DiffractionPatterns.calc_dqshift_for_given_th(angles_list,ki,kf,qbragg);
            
            [X_recip,Y_recip,Z_recip, qz_pixel_size] = DiffractionPatterns.calc_reciprocal_space(dqshift,Npix,depth,d2_bragg);
            
            fft_final = fftn(object);
            fft_shift = fft_final.*exp(-1i*(shift(1)*X_recip + shift(2)*Y_recip + shift(3)*Z_recip));
            
            
            object_final_shift = ifftn(fft_shift);
            
        end
        
        
         function [err_final] = calculate_error_realspace(rho_original,rho_final)
            
            error_3D = (rho_final-rho_original);
            error_3D_mod = error_3D(:)'*error_3D(:);
            err_final = error_3D_mod/numel(error_3D);
            
        end
      
        
    end
end
    

