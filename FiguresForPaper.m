classdef FiguresForPaper
   % This library of function generates the figures of the optic letters
    
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [rho_3DFT] = effect_angle_shift_3DFT(NW,probe,delta_thscanvals,index_to_distort,percent,ki_o,kf_o,X,Y,Z)
           
            [dth_disp] = Phretrieval_functions.generate_angular_jitter(delta_thscanvals,index_to_distort,percent);
            
            %[dq_shift_nominal] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,kf_o-ki_o);
            
            [dq_shift_real] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals + dth_disp',ki_o,kf_o,kf_o-ki_o);
            
            [ simI,rock_curve,Proj_vol,FT_Psij,Qterm] = DiffractionPatterns.calc_dp(dq_shift_real,probe,NW,X,Y,Z);

            for ii = 1:numel(FT_Psij)
               Psi(:,:,ii) = FT_Psij(ii).FT_Psij; 
            end
            
            rho_3DFT = fftshift(ifftn(fftshift(Psi)));
            
        end

        function [rho_2DFT] = effect_angle_shift_2DFT(NW,probe,delta_thscanvals,index_to_distort,percent,ki_o,kf_o,X,Y,Z)
            
            [dth_disp] = Phretrieval_functions.generate_angular_jitter(delta_thscanvals,index_to_distort,percent);
            
            [dq_shift_nominal] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,kf_o-ki_o);
            
            [dq_shift_real] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals + dth_disp',ki_o,kf_o,kf_o-ki_o);
            
            [ simI_real,rock_curve_real,Proj_vol_real,FT_Psij_real,Qterm_real] = DiffractionPatterns.calc_dp(dq_shift_real,probe,NW,X,Y,Z);
            
            [Psi_mod_nominal,~,~,FT_Psij_nominal,Qterm_nominal] = DiffractionPatterns.calc_dp(dq_shift_nominal,probe,NW,X,Y,Z);
            
            
            rho_2DFT = zeros(size(NW));
            
            for jj=1:numel(FT_Psij_real)
                
                Pmrho_dummy = fftshift(ifftn(fftshift(FT_Psij_real(jj).FT_Psij)));
                Pmrho_dummy2 = repmat(Pmrho_dummy,[1 1 numel(FT_Psij_real)])/numel(FT_Psij_real);
                rho_2DFT =  rho_2DFT + conj(Qterm_nominal(jj).Qterm).*Pmrho_dummy2;
                
            end
        end
        
        function [rho_2DFT,mn,error_2DFT,data] = effect_angle_shift_2DFT_noise(NW,probe,delta_thscanvals,dth_disp,mncrate,ki_o,kf_o,X,Y,Z)
            
            
            [dq_shift_nominal] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,kf_o-ki_o);
            
            [dq_shift_real] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals + dth_disp',ki_o,kf_o,kf_o-ki_o);
            
            [ simI_real,rock_curve_real,Proj_vol_real,FT_Psij_real,Qterm_real] = DiffractionPatterns.calc_dp(dq_shift_real,probe,NW,X,Y,Z);
            
            [Psi_mod_nominal,~,~,FT_Psij_nominal,Qterm_nominal] = DiffractionPatterns.calc_dp(dq_shift_nominal,probe,NW,X,Y,Z);
            
            
            rho_2DFT = zeros(size(NW));
            
            middpind = round(numel( simI_real)/2);
            mn = mean( simI_real(middpind).Psi_mod(:));
            
            for jj=1:numel(FT_Psij_real)
                data(jj).I = poisrnd(simI_real(jj).Psi_mod*mncrate/mn);
                FT_real(jj).noise_Psij = sqrt( data(jj).I ).*exp(1i*angle(FT_Psij_real(jj).FT_Psij));
                Pmrho_dummy = fftshift(ifftn(fftshift(FT_real(jj).noise_Psij)));
                Pmrho_dummy2 = repmat(Pmrho_dummy,[1 1 numel(FT_Psij_real)])/numel(FT_Psij_real);
                rho_2DFT =  rho_2DFT + conj(Qterm_nominal(jj).Qterm).*Pmrho_dummy2;
                
            end
            
            error_2DFT = DiffractionPatterns.calc_error_multiangle(probe, rho_2DFT, data,delta_thscanvals,ki_o,kf_o,X,Y,Z);
        end
        
        function [h] = display_effect_angle_noise_dp_2D(dp_struct,info_struct,window,slice,dimension,fig_num)
            
            
            
            h = figure(fig_num);
            clf;
            
            set(gcf,'Name',['Count rate = ' num2str(info_struct(1).cntrate) ' Dimension = ' dimension]);
            
            SubplotNum = numel(dp_struct);
            counter = 0;
            
            for jj = 1:numel(dp_struct)
                
                for kk = 1:numel(dp_struct(jj).dp)
                   dp_matrix3D(:,:,kk) =  dp_struct(jj).dp(kk).I;
                end
                
                 switch dimension
                    
                    case '1'
                        NW1_to_plot = squeeze(dp_matrix3D(slice,window(1):window(2),window(1):window(2)));
                    case '2'
                        NW1_to_plot = squeeze(dp_matrix3D(window(1):window(2),slice,window(1):window(2)));
                    case '3'
                        NW1_to_plot = dp_matrix3D(window(1):window(2),window(1):window(2),slice);
                end
                
                
                subplot(SubplotNum,1,jj);
                imagesc(( NW1_to_plot));
                axis image;
                colorbar;
                title([ ' % = ' num2str(info_struct(jj).percent)])
                %caxis([0.022 0.03]);
                
                
            end
            
        end

        
        function [h] = display_effect_angle_noise_2D(rho_struct,info_struct,support,caxis_int,caxis_phase,window,slice,dimension,fig_num)
            
            h = figure(fig_num);            
            clf;
            
            set(gcf,'Name',['Count rate = ' num2str(info_struct(1).cntrate) ' Dimension = ' dimension]);
            
            SubplotNum = numel(rho_struct);
            counter = 0;
            for jj = 1:numel(rho_struct)
                
               rho = rho_struct(jj).rho.*support;
               
                switch dimension
                    
                    case '1'
                        NW1_to_plot = squeeze(rho(slice,window(1):window(2),window(1):window(2)));
                    case '2'
                        NW1_to_plot = squeeze(rho(window(1):window(2),slice,window(1):window(2)));
                    case '3'
                        NW1_to_plot = rho(window(1):window(2),window(1):window(2),slice);
                end
                
                subplot(SubplotNum,2,jj+counter);
                imagesc(abs( NW1_to_plot));
                axis image;
                colorbar;
                title([ ' % = ' num2str(info_struct(jj).percent)])
                caxis([caxis_int]);
                
                subplot(SubplotNum,2,2*jj);
                imagesc(angle( NW1_to_plot));
                axis image;
                colorbar;
                caxis([caxis_phase]);
                
                counter = counter + 1;
            end
            
        end
        
        
        
        function [h] = display_effect_angle_noise_1D(rho_ini,rho_struct,info_struct,window,slice,dimension,fig_num)
            
            h = figure(fig_num);            
            clf;
            
            set(gcf,'Name',['Count rate = ' num2str(info_struct(1).cntrate) ' Dimension = ' dimension]);
            
            switch dimension
                
                case '12'
                    NW1_to_plot = squeeze(rho_ini(slice(1),slice(2),window(1):window(2)));
                case '13'
                    NW1_to_plot = squeeze(rho_ini(slice(1),window(1):window(2),slice(2)));
                case '23'
                    NW1_to_plot = squeeze(rho_ini(window(1):window(2),slice(1),slice(2)));
            end
            
            subplot(1,2,1);
            plot(abs( NW1_to_plot));
            legend_cell{1} = [ 'Original object ' ];
            
            subplot(1,2,2);
            plot(angle( NW1_to_plot));
            
            
            for jj = 1:numel(rho_struct)
                
                switch dimension
                    
                    case '12'
                        NW1_to_plot = squeeze(rho_struct(jj).rho(slice(1),slice(2),window(1):window(2)));
                    case '13'
                        NW1_to_plot = squeeze(rho_struct(jj).rho(slice(1),window(1):window(2),slice(2)));
                    case '23'
                        NW1_to_plot = squeeze(rho_struct(jj).rho(window(1):window(2),slice(1),slice(2)));
                end
                
                subplot(1,2,1);
                hold on;
                plot(abs( NW1_to_plot));
                legend_cell{jj+1} = [ ' % = ' num2str(info_struct(jj).percent)];
                
                subplot(1,2,2);
                hold on;
                plot(angle( NW1_to_plot));
               
                
            end
            subplot(121);
            legend(legend_cell);
            
             subplot(122);
            legend(legend_cell);
            
        end
        
        function [rho_3DFT,FT_set_slice,FTmod_noise,mncrate_complete] = calculate_rho_dp_Figure1(delta_thscanvals,angshift_percentage,mncrate,rho,probe,ki,kf,X,Y,Z)
            % This function does the b) panel of Figure 1 where the detrimental effects
            % of angular uncertainties are shown as a function of the
            % mncrate and the angular shift
            
            % initialize
            dth_disp = zeros(numel(delta_thscanvals),1);
            index_to_distort = [34:numel(delta_thscanvals)];
            FT_set =zeros(size(X));
            FTmod =zeros(size(X));
            support = abs(rho);           
                        
           
            
            for ll = 1:numel(angshift_percentage)
                % calculate the angular shift
                dth_disp(index_to_distort) = -[angshift_percentage(ll)*unique(diff(delta_thscanvals))/100];
                angles_list = delta_thscanvals + dth_disp';
                
                % calculate a list of momentum transfers:
                [dqshift] = DiffractionPatterns.calc_dqshift_for_given_th(angles_list,ki,kf,kf-ki);
                
                % calculate the corresponding set of diffracted waves:
                [Psi_mod,~,~,FT_Psij,~] = DiffractionPatterns.calc_dp(dqshift,probe,rho,X,Y,Z);
                
                % calculate the mean of the medium slice
                middpind = round(numel(Psi_mod)/2);
                mn = mean(Psi_mod(middpind).Psi_mod(:));
              
                mncrate_complete = [mncrate mn];
                for kk = 1:numel(mncrate)
                    for jj =1:numel(FT_Psij)
                        FTmod(:,:,jj) = poisrnd(Psi_mod(jj).Psi_mod*mncrate_complete(kk)/mn);
                        FT_set(:,:,jj) = sqrt(FTmod(:,:,jj)).*exp(1i*angle(FT_Psij(jj).FT_Psij));
                    end
                    rho_3DFT.angshift(ll).mncrate(kk).matrix = fftshift(ifftn(fftshift(FT_set)))*sqrt(mn/mncrate_complete(kk)).*support;
                    FT_set_slice.angshift(ll).mncrate(kk).matrix =  FT_set;
                    FTmod_noise.angshift(ll).mncrate(kk).matrix = FTmod;
                    
                    %savestring = ['NW_perfect_noise_mcnrate' num2str(kk) '_angshift' num2str(ll)];
                    %save(savestring,'FTmod','FT_set');
                end
                
                    for jj =1:numel(FT_Psij)
                        FTmod(:,:,jj) = Psi_mod(jj).Psi_mod;
                        FT_set(:,:,jj) = sqrt(FTmod(:,:,jj)).*exp(1i*angle(FT_Psij(jj).FT_Psij));
                    end
                    rho_3DFT.angshift(ll).mncrate(kk+1).matrix = fftshift(ifftn(fftshift(FT_set)));
                    FT_set_slice.angshift(ll).mncrate(kk+1).matrix =  FT_set;
                    FTmod_noise.angshift(ll).mncrate(kk+1).matrix = FTmod;
                
            end
            
            
        end
        
         
        function [] = display_Figure1(matrix_3DFT,angshift_percentage,mncrate,index_to_plot,fig_num)
            
            
            
            figure(fig_num);
            clf;
            

            for ll = 1:numel(angshift_percentage)
                for kk = 1:numel(matrix_3DFT.angshift(ll).mncrate)
                    
                    %load([loadstring num2str(kk) '_angshift' num2str(ll) '.mat']);
                    
                    middpind = index_to_plot;
                    
                    subplot(numel(angshift_percentage),numel(matrix_3DFT.angshift(ll).mncrate),numel(matrix_3DFT.angshift(ll).mncrate)*(ll-1)+kk);
                    imagesc(abs(matrix_3DFT.angshift(ll).mncrate(kk).matrix(:,:,middpind)));
                    axis image;
                    colorbar;
                    if kk~=3
                        title(['Count rate = ' num2str(mncrate(kk)) ' Angular shift = ' num2str(angshift_percentage(ll))]);
                    else
                        title(['Original object;Angular shift = ' num2str(angshift_percentage(ll))]);                        
                    end
                    
                end
                                
            end
            
            
             figure(fig_num+1);
            clf;
            
             for ll = 1:numel(angshift_percentage)
                for kk = 1:numel(matrix_3DFT.angshift(ll).mncrate)
                    
                    %load([loadstring num2str(kk) '_angshift' num2str(ll) '.mat']);
                    
                    middpind = index_to_plot;
                    
                    subplot(numel(angshift_percentage),numel(matrix_3DFT.angshift(ll).mncrate),numel(matrix_3DFT.angshift(ll).mncrate)*(ll-1)+kk);
                    imagesc(angle(matrix_3DFT.angshift(ll).mncrate(kk).matrix(:,:,middpind)));
                    axis image;
                    colorbar;
                     if kk~=3
                        title(['Count rate = ' num2str(mncrate(kk)) ' Angular shift = ' num2str(angshift_percentage(ll))]);
                    else
                        title(['Original object;Angular shift = ' num2str(angshift_percentage(ll))]);                        
                    end
                    
                end
                                
            end
            
            
        end
        
        function [Psi_mod_matrix] = display_slice_dp(rho,probe,angles_list,slice_array ,ki,kf,d2_bragg,thBragg,X,Y,Z,fig_num)
            
            Npix = size(X,1);
            depth = size(X,3);
            
            qbragg = kf-ki;
            % calculate a list of momentum transfers:
            [dqshift,ki_list,kf_list] = DiffractionPatterns.calc_dqshift_for_given_th(angles_list,ki,kf,qbragg);
            
            % pixel size on the z direction of the reciprocal space:
            qz_pixel_size = abs((dqshift(1,3) - dqshift(end,3))/size(dqshift,1));
            
            % calculate the corresponding set of diffracted waves:
            [Psi_mod,~,~,FT_Psij,~] = DiffractionPatterns.calc_dp(dqshift,probe,rho,X,Y,Z);
            
            for jj = 1:numel(Psi_mod)
               Psi_mod_matrix(:,:,jj) = Psi_mod(jj).I; 
            end
            
            figure(fig_num);
            clf;
            
            subplot(121);
            hold on;
            quiver3(0,0,0,ki(1),ki(2),ki(3),'k')
            quiver3(0,0,0,kf(1),kf(2),kf(3),'k')
            quiver3(0,0,0,kf(1)-ki(1),kf(2)-ki(2),kf(3)-ki(3),'k')
            
            jj = 1;
            
            Rvect = [0 1 0];
            Rcenter = [0 0 0];
            
            % plane 1 of the beam
            [ grid(jj).yy1,  grid(jj).xx1 ,  grid(jj).zz1] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), [0]);
            
            % sets the cartesian coordinates
            grid(jj).points2d_1 = [  grid(jj).xx1(:)  grid(jj).yy1(:)];
            
            
            % create a 3rd dimension and convert the pixel number in inverse angstroms:
            grid(jj).points3d_1 = [  grid(jj).points2d_1 grid(jj).zz1(:)];
            
            
            % reshape the mesh grid which is now centered around the beam maximum
            % intensity and in microns units:
            grid(jj).xx1_microns = reshape(grid(jj).points3d_1(:,1),size(grid(jj).xx1));
            grid(jj).yy1_microns = reshape(grid(jj).points3d_1(:,2),size(grid(jj).xx1));
            grid(jj).zz1_microns = reshape(grid(jj).points3d_1(:,3),size(grid(jj).xx1));
            
            % plot
            subplot(122);
            hold on;
            h(jj) = surf(grid(jj).xx1_microns, grid(jj).yy1_microns, grid(jj).zz1_microns);
            h(jj).EdgeColor = 'none'; % set properties of the plot
            h(jj).CData = abs(Psi_mod_matrix(:,:,jj)./max(Psi_mod_matrix(:)));
            alpha(h(jj),0.3);
                
            axis image
            
            % rotation
            rotate(h(jj), Rvect, thBragg, Rcenter);
            
            hold off;
            
            % extract the rotated coordinates in the lab frame:
            
            rotpos.x = h(jj).XData;
            rotpos.y = h(jj).YData;
            rotpos.z = h(jj).ZData;
            
            
            % fft space:
            [X_recip,Y_recip,Z_recip] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg),...
                [-depth/2:depth/2-1].*2*qz_pixel_size);
            
            
            hold on;
            h_Psi2 =  di(Psi_mod_matrix./max(Psi_mod_matrix(:)),2e-4,'g',X_recip,Y_recip,Z_recip);
            
%             [X_recip_toplot Y_recip_toplot Z_recip_toplot] = meshgrid([-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg), ...
%                 [-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg),...
%                 [-depth/2 depth/2-1].*2*pi/(depth*d2_bragg));
%             
%             scatter3(X_recip_toplot(:),Y_recip_toplot(:),Z_recip_toplot(:));
            
            
            for jj = slice_array 
                
                subplot(121);
                hold on;
                quiver3(0,0,0,ki_list(jj,1),ki_list(jj,2),ki_list(jj,3),'r');
                quiver3(0,0,0,kf_list(jj,1),kf_list(jj,2),kf_list(jj,3),'b');
                quiver3(0,0,0,kf_list(jj,1)-ki_list(jj,1),kf_list(jj,2)-ki_list(jj,2),kf_list(jj,3)-ki_list(jj,3),'g')
                quiver3(qbragg(1),qbragg(2),qbragg(3),dqshift(jj,1),dqshift(jj,2),dqshift(jj,3),'r')
                
                %Psi2_toplot(:,:,jj) = temp;
                %translate the detector plane from dq_shift:
                trrotpos(jj).x = rotpos.x + dqshift(jj,1);
                trrotpos(jj).y = rotpos.y + dqshift(jj,2);
                trrotpos(jj).z = rotpos.z + dqshift(jj,3);
                
                
                % plot
                subplot(122);
                hold on;
                h(jj) = surf(trrotpos(jj).x, trrotpos(jj).y, trrotpos(jj).z);
                h(jj).EdgeColor = 'none'; % set properties of the plot
                h(jj).CData = abs(Psi_mod_matrix(:,:,jj)./max(Psi_mod_matrix(:))); % use the intensity of the beam to color the surface
                alpha(h(jj),0.3);                
                drawnow;
                
                axis image;
                %hold off;
                
                %quiver3(rotpos.x(1),rotpos.y(1),rotpos.z(1),dqshift(jj,1),dqshift(jj,2),dqshift(jj,3), 'r','LineWidth',3,'ShowArrowHead','on');

            end
            
            
        end
        
        function [] = figure2_rightpanel(NW1,NW2,NW3,titleNW1,titleNW2,titleNW3,intenscolor,phasecolor,window,index,dimension,fig_num)
            
            figure(fig_num);
            clf;
            
            
            
            for ii = index
                
                switch dimension
                    
                    case '1'
                        NW1_to_plot = squeeze(NW1(ii,window(1):window(2),window(3):window(4)));
                        NW2_to_plot = squeeze(NW2(ii,window(1):window(2),window(3):window(4)));
                        NW3_to_plot = squeeze(NW3(ii,window(1):window(2),window(3):window(4)));
                    case '2'
                        NW1_to_plot = squeeze(NW1(window(1):window(2),ii,window(3):window(4)));
                        NW2_to_plot = squeeze(NW2(window(1):window(2),ii,window(3):window(4)));
                        NW3_to_plot = squeeze(NW3(window(1):window(2),ii,window(3):window(4)));
                    case '3'
                        NW1_to_plot = NW1(window(1):window(2),window(3):window(4),ii);
                        NW2_to_plot = NW2(window(1):window(2),window(3):window(4),ii);
                        NW3_to_plot = NW3(window(1):window(2),window(3):window(4),ii);
                end
                
                
                subplot(321);
                imagesc(abs(NW1_to_plot));
                axis image;
                colorbar;
                title([titleNW1]);
                set(gca,'FontSize',30);
                axis off;
                caxis(intenscolor);
                
                subplot(322);
                imagesc(angle(NW1_to_plot));
                axis image;
                colorbar;
                title([titleNW1 ]);
                set(gca,'FontSize',30);
                axis off;
                caxis(phasecolor);
                
                subplot(323);
                imagesc(abs(NW2_to_plot));
                axis image;
                colorbar;
                title([titleNW2 ]);
                set(gca,'FontSize',30);
                axis off;
                caxis(intenscolor);
                
                subplot(324);
                imagesc(angle(NW2_to_plot));
                axis image;
                colorbar;
                title([titleNW2 ]);
                set(gca,'FontSize',30);
                axis off;
                caxis(phasecolor);
                
                
                subplot(325);
                imagesc(abs(NW3_to_plot));
                axis image;
                colorbar;
                title([titleNW3 ]);
                set(gca,'FontSize',30);
                axis off;
                caxis(intenscolor);
                
                subplot(326);
                imagesc(angle(NW3_to_plot));
                axis image;
                colorbar;
                title([titleNW3 ]);
                set(gca,'FontSize',30);
                axis off;
                caxis(phasecolor);
                
                pause(.5);
                
                
                
            end
            
            
            
            
        end
        
        
         function [] = figure5_bottompanel(struct_toplot,intenscolor,phasecolor,phasecolor_2,window,index,dimension,fig_num)
            
            figure(fig_num);
            clf;
            
            
            
            for jj = 1:numel(struct_toplot)
                
                NW1_3D = struct_toplot(jj).rho_shift;
                support_3D = struct_toplot(jj).support_shift_fin;
                
                NW1_3D_nophase = struct_toplot(jj).rho_nophase;
                
                for ii = index
                    
                    switch dimension
                        
                        case '1'
                            NW1_to_plot = squeeze(NW1_3D(ii,window(1):window(2),window(3):window(4)));
                            NW1_nophase = squeeze(NW1_3D_nophase(ii,window(1):window(2),window(3):window(4)));
                        case '2'
                            NW1_to_plot = squeeze(NW1_3D(window(1):window(2),ii,window(3):window(4)));
                            NW1_nophase = squeeze(NW1_3D_nophase(window(1):window(2),ii,window(3):window(4)));
                        case '3'
                            NW1_to_plot = NW1_3D(window(1):window(2),window(3):window(4),ii);
                            NW1_nophase = NW1_3D_nophase(window(1):window(2),window(3):window(4),ii);
                    end
                    

                    
                    
                    subplot(3,numel(struct_toplot),jj);
                    imagesc(abs(NW1_to_plot));
                    axis image;
                    colorbar;          
                    set(gca,'FontSize',30);
                    axis off;
                    caxis(intenscolor);
                    
                    subplot(3,numel(struct_toplot),numel(struct_toplot)+jj);
                    imagesc(angle(NW1_to_plot));
                    axis image;
                    colorbar;
                    set(gca,'FontSize',30);
                    axis off;
                    caxis(phasecolor);
                    
                    subplot(3,numel(struct_toplot),2*numel(struct_toplot)+jj);
                    imagesc(angle(NW1_nophase));
                    axis image;
                    colorbar;
                    set(gca,'FontSize',30);
                    axis off;
                    caxis(phasecolor_2);
                    
                    
                    pause(.5);
                    
                    
                    
                end
                
            end
            
            
        end
        
        
        
        function [support_new_shift_final,support_new_shift_vector] = figure2_flipsupport(flipflag,NW,support_new,angle_list,probe,ki_o,kf_o,X,Y,Z,d2_bragg)
            
            support_new_2DFT_ini = DiffractionPatterns.From3DFT_to_2DFT(support_new,angle_list,probe,ki_o,kf_o,X,Y,Z);
            
            
            if flipflag == 1
                support_new_2DFT = ifftn(conj(fftn(support_new_2DFT_ini)));
            else
                support_new_2DFT = support_new_2DFT_ini;
            end
            
            [support_new_shift,support_new_shift_vector] = DiffractionPatterns.shift_object(NW,support_new_2DFT,angle_list,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
            
            support_shift_abs = abs(support_new_shift);
            support_new_shift_final = (support_shift_abs>0.1*max(support_shift_abs(:)));
            
            
            
            DisplayResults.compare_two_objects(NW,support_new_shift_final,'','',[40 90 40 90],[65],'3',31);
            
        end
        
    end
    
    
end