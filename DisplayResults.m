classdef DisplayResults
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [theta_iter] = read_angles_iterations(data,ini_angle_val,delta_thscanvals)
            
            % read values for the different angles and the iterations
            for jj = 1:numel(data)
                
                % iteration 0
                kk = 1;
                
                theta_iter(jj).dth_new_iter(kk) = data(jj).theta_iter_ini.dth_new_iter;
                
                difftheta = unique(diff(delta_thscanvals));
                
                if jj>1
                    theta_iter(jj).dth_new_neighbor(kk) = (data(jj).theta_iter_ini.dth_new_iter -  data(jj-1).theta_iter_ini.dth_new_iter)*100/difftheta(1);
                else
                    theta_iter(jj).dth_new_neighbor(kk) = (data(jj+1).theta_iter_ini.dth_new_iter - data(jj).theta_iter_ini.dth_new_iter)*100/difftheta(1);
                end
                
                theta_iter(jj).beta(kk) = data(jj).theta_iter_ini.beta;
                theta_iter(jj).grad_final(kk) = data(jj).theta_iter_ini.grad_final;
                theta_iter(jj).dqshift(kk,:) = data(jj).theta_iter_ini.dqshift;
                theta_iter(jj).delta_dth(kk) = data(jj).theta_iter_ini.dth_new_iter-ini_angle_val(jj);
                
                
                
                
                for kk = 1:numel(data(jj).theta_iter)
                    theta_iter(jj).dth_new_iter(kk+1) = data(jj).theta_iter(kk).dth_new_iter;
                    
                    if jj>1
                        theta_iter(jj).dth_new_neighbor(kk+1) = (data(jj).theta_iter(kk).dth_new_iter -  data(jj-1).theta_iter(kk).dth_new_iter)*100/difftheta(1);
                    else
                        theta_iter(jj).dth_new_neighbor(kk+1) = (data(jj+1).theta_iter(kk).dth_new_iter - data(jj).theta_iter(kk).dth_new_iter)*100/difftheta(1);
                    end
                    
                    theta_iter(jj).beta(kk+1) = data(jj).theta_iter(kk).beta;
                    theta_iter(jj).grad_final(kk+1) = data(jj).theta_iter(kk).grad_final;
                    theta_iter(jj).dqshift(kk+1,:) = data(jj).theta_iter(kk).dqshift;
                    theta_iter(jj).delta_dth(kk+1) = data(jj).theta_iter(kk).dth_new_iter-ini_angle_val(jj);
                end
            end
            
            
        end
               
        function [h] = display_angles_iterations(theta_iter,data,angles_to_plot,niter,fig_ini)
            
            fig_num = [1:numel(angles_to_plot)]+fig_ini;
            
            for ff = 1:numel(fig_num)
                figure(fig_num(ff))
                clf;
                
                for jj = 1:numel(theta_iter)
                    
                    if jj == angles_to_plot(ff)%mod(jj,5) == 0
                        subplot(3,2,1);
                        hold on;
                        plot(theta_iter(jj).dth_new_iter(1:niter),'-.');
                        plot(data(jj).dth_real*ones(niter,1),'r');
                        xlabel('iterations');ylabel('theta');
                        title(['theta vs iterations @ theta = ' num2str(data(jj).dth_nominal)]);
%                         
%                         hold on;
%                         plot((theta_iter(jj).dth_new_iter(1:niter)-data(jj).dth_real)./data(jj).dth_real,'-.');
%                         xlabel('iterations');ylabel('theta');
%                         title(['relative error @ theta = ' num2str(data(jj).dth_nominal)]);
%                         
%                         
                        subplot(3,2,2);
                        hold on;
                        plot(theta_iter(jj).beta(1:niter))
                        xlabel('iterations');ylabel('beta');
                        title('beta vs iterations');
                        
                        subplot(3,2,3);
                        hold on;
                        plot(theta_iter(jj).grad_final(1:niter))
                        xlabel('iterations');ylabel('gradient');
                        title('gradient vs iterations');
                        
                        
                        subplot(3,2,4);
                        hold on;
                        plot(theta_iter(jj).dth_new_neighbor(1:niter));
                        plot(theta_iter(jj).dth_new_neighbor(1)*ones(numel(data(jj).theta_iter),1),'r');

                        xlabel('iterations');ylabel('\Delta\theta neighbors');
                        title('Percentage of displacement with respect to neighbor');
                        
                        
                        subplot(3,2,5);
                        hold on;
                        plot(theta_iter(jj).delta_dth(1:niter))
                        xlabel('iterations');ylabel('delta_dth');
                        title('delta_dth vs iterations');
                        
                        
                        subplot(3,2,6);
                        hold on;
                        plot(theta_iter(jj).dqshift(1:niter,1),'.');
                        plot(data(jj).dqshift_real(1)*ones(numel(data(jj).theta_iter(1:niter)),1),'r');
                      
                        xlabel('iterations');ylabel('dq components');
                        title('dqshift vs iterations');
                        legend('qx','real qx');%,'qy','real qy','qz','real qz');
                        
                        set(gcf,'Name',['th = ' num2str(data(jj).dth_nominal)]);

                        
                    end
                    
                end
                
                
            end
            
            
        end
       
        function [h] = display_angles_iterations_errorrel(theta_iter,data,dth_disp,angles_to_plot,niter,flag,fig_ini)
            
            fig_num = [1:numel(angles_to_plot)]+fig_ini;
            NumofSubplots = numel(angles_to_plot);
            
            figure(fig_ini);
            clf;
            
            for ff = 1:numel(fig_num)
                
                
                for jj = 1:numel(theta_iter)
                                                            
                    if jj == angles_to_plot(ff)%mod(jj,5) == 0
                        subplot(ceil(sqrt(NumofSubplots)),ceil(sqrt(NumofSubplots)),ff);
                        hold on;
                        switch flag
                            case 'relative'
                                %plot((theta_iter(jj).dth_new_iter(1:niter)-data(jj).dth_real)./data(jj).dth_real,'-.');
                                error_rel = (theta_iter(jj).dth_new_iter(1:niter)-data(jj).dth_nominal -dth_disp(jj))./dth_disp(jj);
                                plot(error_rel,'.-','LineWidth',3.0);
                                hold on;
                                plot(zeros(niter,1),'r');
                                ylabel('Relative error ');
                                Namestr = ['Relative error'];
                            case 'absolute'
                                plot((theta_iter(jj).dth_new_iter(1:niter)),'-.','LineWidth',3.0);
                                hold on;
                                plot(data(jj).dth_real*ones(niter,1),'r');
                                ylabel('angle');
                                 Namestr = ['Absolute error'];
                        end
                        xlabel('iterations');
                        title(['th = ' num2str(data(jj).dth_nominal)]);
                        set(gcf,'Name',Namestr);

                        
                    end
                    
                end
                
                
            end
            
            
         end
        
        function [h] = display_all_angles_oneiterations_errorrel(theta_iter,data,dth_disp,nitervect,flag,fig_ini)
             
             figure(fig_ini);
             clf;
             hold on;
            for kk = 1:numel(nitervect)
             for jj = 1:numel(data)
                 dth_nom(jj) = data(jj).dth_nominal;
                 hold on;
                 switch flag
                     case 'relative'
                         %plot((theta_iter(jj).dth_new_iter(1:niter)-data(jj).dth_real)./data(jj).dth_real,'-.');
                         
                         error_rel(kk).err(jj) = (theta_iter(jj).dth_new_iter(nitervect(kk))-data(jj).dth_nominal -dth_disp(jj))./dth_disp(jj);
                         Namestr = ['Relative error'];

                     case 'absolute'
                         error_rel(kk).err(jj) = theta_iter(jj).dth_new_iter(nitervect(kk))-data(jj).dth_nominal -dth_disp(jj);
                         hold on;
                         %plot(data(jj).dth_real*zeros(numel(data),1),'r');
                         ylabel('angle');
                         Namestr = ['Absolute error'];

                 end
             end
             
             plot(dth_nom,error_rel(kk).err,'.-','MarkerSize',15);
             
             
             legend_str{kk} = ['iter = ' num2str(nitervect(kk))];

             set(gcf,'Name', Namestr);
            end
            legend(legend_str);
            plot(dth_nom,zeros(numel(data),1),'r','LineWidth',5);
             ylabel('Angular error in degrees');
            xlabel('Nominal angle in degrees');
            set(gca,'FontSize',30);
            set(gca,'LineWidth',2);
           

         end
        
        
        function [h] = display_neighboring_angles_iterations(theta_iter,data,angles_to_plot,flag,fig_ini)
            
            fig_num = [1:numel(angles_to_plot)]+fig_ini;
            NumofSubplots = numel(angles_to_plot);
            
            figure(fig_ini);
            clf;
            
            for ff = 1:numel(fig_num)
                
                
                for jj = 1:numel(theta_iter)
                                                            
                    if jj == angles_to_plot(ff)%mod(jj,5) == 0
                        subplot(ceil(sqrt(NumofSubplots)),ceil(sqrt(NumofSubplots)),ff);
                        hold on;
                        
                        switch flag
                            case 'relative'
                                plot(theta_iter(jj).dth_new_neighbor - theta_iter(jj).dth_new_neighbor(1));
                            case 'absolute'
                                plot(theta_iter(jj).dth_new_neighbor);                             
                                plot(theta_iter(jj).dth_new_neighbor(1)*ones(numel(data(jj).theta_iter),1),'r');
                        end
                        xlabel('iterations');ylabel('%\Delta \theta\_j+1 - \theta\_j ');
                        title(['th = ' num2str(data(jj).dth_nominal)]);
                        set(gcf,'Name',['Percentage of displacement with respect to neighbor']);

                        
                    end
                    
                end
                
                
            end
            
            
        end
        
        function [] = show_rho_theta_update(figure_num,errlist,rho,midsl,anglelist,trueangles,norm_grad_rho,beta_rho,norm_grad_theta,beta_theta,flag)
            % this function creates the figure where we plot the improvements of the
            % rho and the angles during the phase retrieval iterations. flagINI
            % indicates if the figure needs to be created (initial state) or only
            % updated
            midsl_2 = round(size(rho,1)/2);
            
            figure(figure_num);
            if strcmp(flag,'Ini')
                clf;
                setfigsize(gcf, 1000,500);
            end
            % plot
            subplot(241); 
            %imagecomp(rho(:,:,midsl)); 
            imagesc(abs(rho(:,:,midsl))); 
            colorbar; axis image;title('Intensity'); %zoom(1.5);
            
            subplot(242); 
            imagesc(angle(rho(:,:,midsl))); 
            %imagecomp(squeeze(rho(100,:,:))); 
            colorbar; axis image;title('Phase') %zoom(1.5);
            
            subplot(243); plot(log10(errlist),'LineWidth',3.0);title('Error metric');xlabel('Iterations')
            
            
            subplot(244);cla;
            plot(trueangles,(anglelist-trueangles)./trueangles,'ob');
            hold on;
            plot(trueangles,zeros(numel(trueangles),1),'r');
            title('Relative error in Angles');
            xlabel('angle (deg)');
           % hold on; plot(trueangles,'*r');
            
               
              subplot(245); 
            %imagecomp(rho(:,:,midsl)); 
            imagesc(abs(squeeze(rho(midsl_2,:,:)))); 
            colorbar; axis image;title('Intensity'); %zoom(1.5);
            subplot(246); 
            imagesc(angle(squeeze(rho(midsl_2,:,:)))); 
            %imagecomp(squeeze(rho(100,:,:))); 
            colorbar; axis image;title('Phase') %zoom(1.5);
            
            
            subplot(247);
            yyaxis left;
            plot(log10(norm_grad_rho),'ob');
            title(['Rho, grad(end) = ' num2str(norm_grad_rho(end))]);
            xlabel('iterations');
            
            ax = gca;
            set(ax.YAxis(1),'Color','b');
            ylabel('log norm of gradien');
            
            yyaxis right;
            ax = gca;
            set(ax.YAxis(2),'Color','k');
            plot(beta_rho,'*k');
            ylabel('beta');
           
            subplot(248);
            yyaxis left;
            plot(log10(norm_grad_theta),'ob');
            if numel(norm_grad_theta) == 0
                title(['Theta grad(end) = ' num2str(norm_grad_theta)]);
            else
                title(['Theta grad(end) = ' num2str(norm_grad_theta(numel(norm_grad_theta)))]);                
            end
            
            xlabel('iterations');
            ylabel('log norm of gradien');
            
            ax = gca;
            set(ax.YAxis(1),'Color','b');
            
            yyaxis right;
            ax = gca;
            set(ax.YAxis(2),'Color','k');
            plot(beta_theta,'*k');
            ylabel('beta');
            drawnow;
            
        end
        
        function [h,lgd] = show_rocking_curve(angles,rock_curve,flag,fig_num,color,legend_str)
           % this function displays the rocking curve with respect to a specified angular grid
           
           if strcmp(flag,'new')
               figure(fig_num);
               clf;
               h = bar(angles,sqrt(rock_curve),'FaceColor',color,'BarWidth',0.2,'EdgeColor',color);
               lgd = legend(legend_str);
           elseif strcmp(flag,'hold')
               figure(fig_num);
               ax = gca;
               hold on;
               h = bar(angles,sqrt(rock_curve),'FaceColor',color,'BarWidth',0.2,'EdgeColor',color);
               ax.Legend.String{size(ax.Legend.String,2)} = legend_str;
               lgd = legend(ax.Legend.String);
           end
           
           xlabel('Rock angle');
           ylabel('sqrt(Psi.*conj(Psi))');
           set(gca,'FontSize',20.0);
           
        end
              
        
        function [] = compare_dpcalc_vs_dpsim(data_exp,data_calc,index,delta_thscanvals,figure_num)
            
            figure(figure_num);
            clf; setfigsize(gcf, 800,400);
            colormap jetvar;
            
            ca = [0 2];
            for ii=index
                
                subplot(141);
                imagesc( sqrt(data_exp(ii).simI)); axis image;colorbar;
                title(['simulated data no noise ii = ' num2str(ii) ' , dth = ' num2str(delta_thscanvals(ii))]);
                
                
                subplot(142);
                imagesc( sqrt(data_exp(ii).I));axis image;colorbar;
                title('simulated data + noise')
                
               
                subplot(143)
                imagesc( sqrt(data_calc(ii).Psi_mod));axis image;colorbar;
                title('calculated data from retrieved object');
                
                 subplot(144)
                imagesc( (sqrt((data_exp(ii).I)) - sqrt(data_calc(ii).Psi_mod)).^2);axis image;colorbar;
                title('difference squared');
                
                 drawnow;
                
                pause();
            end
            
        end
        
        function [] = compare_2DFTdp_vs_3DFTdp(data_3D,data_2D,index,delta_thscanvals,figure_num)
            
            figure(figure_num);
            clf; setfigsize(gcf, 800,400);
            colormap jetvar;
            
            ca = [0 2];
            for ii=index
                
                subplot(131);
                imagesc( sqrt(data_3D(:,:,ii))); axis image;colorbar;
                title(['simulated data no noise ii = ' num2str(ii) ' , dth = ' num2str(delta_thscanvals(ii))]);
                
                                            
               
                subplot(132)
                imagesc( sqrt(data_2D(ii).Psi_mod));axis image;colorbar;
                title('calculated data from retrieved object');
                
                 subplot(133)
                imagesc( (sqrt((data_3D(:,:,ii))) - sqrt(data_2D(ii).Psi_mod)).^2);axis image;colorbar;
                title('difference squared');
                
                 drawnow;
                
                pause(.5);
            end
            
        end
        
        
        function [] = compare_two_objects(NW1,NW2,titleNW1,titleNW2,window,index,dimension,fig_num)
            
            figure(fig_num);
            clf;
            
            
            if numel(dimension) == 1
                for ii = index
                    
                    switch dimension
                        
                        case '1'
                            NW1_to_plot = squeeze(NW1(ii,window(1):window(2),window(3):window(4)));
                            NW2_to_plot = squeeze(NW2(ii,window(1):window(2),window(3):window(4)));
                        case '2'
                            NW1_to_plot = squeeze(NW1(window(1):window(2),ii,window(3):window(4)));
                            NW2_to_plot = squeeze(NW2(window(1):window(2),ii,window(3):window(4)));
                        case '3'
                            NW1_to_plot = NW1(window(1):window(2),window(3):window(4),ii);
                            NW2_to_plot =NW2(window(1):window(2),window(3):window(4),ii);
                    end
                    
                    
                    subplot(221);
                    imagesc(abs(NW1_to_plot));
                    axis image;
                    colorbar;
                    title([titleNW1]);
                    set(gca,'FontSize',30);
                    axis off;
                    
                    subplot(222);
                    imagesc(angle(NW1_to_plot));
                    axis image;
                    colorbar;
                    title([titleNW1 ]);
                    set(gca,'FontSize',30);
                     axis off;
                     
                    subplot(223);
                    imagesc(abs(NW2_to_plot));
                    axis image;
                    colorbar;
                    title([titleNW2 ]);
                    set(gca,'FontSize',30);
                    axis off;
                    
                    subplot(224);
                    imagesc(angle(NW2_to_plot));
                    axis image;
                    colorbar;
                    title([titleNW2 ]);
                    set(gca,'FontSize',30);
                    axis off;
                    
                    pause(.5);
                    
                    
                    
                end
                
            else
                
                ii = index;
                
                switch dimension
                    
                    case '12'
                        NW1_to_plot = squeeze(NW1(ii(1),ii(2),:));
                        NW2_to_plot = squeeze(NW2(ii(1),ii(2),:));
                    case '13'
                        NW1_to_plot = squeeze(NW1(ii(1),:,ii(2)));
                        NW2_to_plot = squeeze(NW2(ii(1),:,ii(2)));
                    case '23'
                        NW1_to_plot = squeeze(NW1(:,ii(1),ii(2)));
                        NW2_to_plot = squeeze(NW2(:,ii(1),ii(2)));
                end
                
                figure(fig_num);
                clf;
                
                subplot(121);
                plot(abs(NW1_to_plot));
                hold on;
                plot(abs(NW2_to_plot),'LineStyle','-.');
                legend(titleNW1,titleNW2);
                title(['Itensity dimensions = ' dimension]);
                set(gca,'FontSize',30);
                
                 subplot(122);
                plot(angle(NW1_to_plot));
                hold on;
                plot(angle(NW2_to_plot));
                legend(titleNW1,titleNW2);
                title(['Phase dimensions = ' dimension]);
                set(gca,'FontSize',30);
                
            end
            
            
        end
        
    end
end