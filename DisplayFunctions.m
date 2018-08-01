classdef DisplayFunctions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function h = display_diff_geom(rho,ki,kf,qbragg,fig_num,X,Y,Z)
            
            if ishandle(fig_num)
                handleFLAG = get(fig_num,'type');
                switch handleFLAG
                    case 'figure'
                        figure(fig_num);
                        
                    case 'axes'
                        subplot(fig_num);
                end
            else
                figure(fig_num)
            end
            
            clf;
            
            hold on;
            
            % the sample 
            %h=di(angle(rho), -.5, 'y', X,Y,Z); alpha(h,.5);
            isosurface( ...
                X,Y,Z, ...
                smooth3( abs(rho), 'gaussian', 9 ), .9, ...
                smooth3( angle( rho ), 'gaussian', 9 ) );
            
            % the wavevectors
%             quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
%             quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
%             quiver3(0,0,0, qbragg(1), qbragg(2), qbragg(3), 'b');
%             
            % the detector 
%             [Xd Yd] = meshgrid([-.1 .1]);
%             surf(Xd,Yd,ones(size(Xd)));
%             view(-2,53);
            
            xlabel('x');
            ylabel('y'); 
            zlabel('z');
            %axis image;
            axis square; 
            axis tight; 
            axis off; 
            colorbar; 
            lighting gouraud;

            %colormap default
            %colorbar
            title('units of displacement (microns)')
                      
        end
        
        function h = display_strain(strptsmid,phvals,fig_num)
            
             
             if ishandle(fig_num)
                 handleFLAG = get(fig_num,'type');
                 switch handleFLAF
                     case 'figure'
                         figure(fig_num);
                         
                     case 'axes'
                         subplot(fig_num);
                 end
             else
                 figure(fig_num)
             end
             
             
            clf;
            
            hold on;
            
            % the strain
            h = scatter3(strptsmid(:,1), strptsmid(:,2), strptsmid(:,3), [], phvals(:));

                      
        end
        
        function h2 = display_calc_dqshift_deriv(dth_nominal,dq_shift_deriv0,dq_shift0,ki,kf,h1,legend_str)
           
            % h1 is the handle for the figure            
            th_fine_grid  = dth_nominal+[-180:10:180];                                
                                                                                        
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(th_fine_grid,ki,kf,kf-ki);
         
            h2 = figure(h1);
            clf;
            
            subplot(1,2,1);
            hold on;
            plot(th_fine_grid,dq_shift(:,1),'b');
            plot(dth_nominal,dq_shift0(1),'ob');
            cstx = dq_shift0(1)-dq_shift_deriv0(1)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
            plot( th_fine_grid,cstx +dq_shift_deriv0(1).*th_fine_grid,'r');
            
            title(['delta q_x at a nominal angle of ' num2str(dth_nominal)])
            legend('dq\_{shift} calculated with the rotation matrix fine grid ','dq\_{shift} calculated with matrix at th_nominal',['derivative of dq\_{shift} calculated ' legend_str])
            
            subplot(1,2,2);
            hold on;
            plot(th_fine_grid,dq_shift(:,3),'b');
            plot(dth_nominal,dq_shift0(3),'ob');
            cstz = dq_shift0(3)-dq_shift_deriv0(3)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
            plot( th_fine_grid, cstz +dq_shift_deriv0(3).*th_fine_grid,'r');
            title(['delta q_z' ]);
           % legend('dq\_{shift} calculated with the rotation matrix ','dq\_{shift} calculated analytically ','derivative of dq\_{shift} calculated analytically')
           
            
        end
 
        function display_calc_grad_theta(probe, rho, data, dth_nominal,grad_calc,err_0 ,X,Y,Z,ki_o,kf_o,fig_num)        
            
            th_fine_grid = dth_nominal+[-0.05:0.5e-3 :0.05];%[-0.005:1e-5 :0.005];
            
            err_thin_grid = zeros(numel(th_fine_grid) ,1);
            
            for jj = 1:numel(th_fine_grid)            
                err_thin_grid(jj) =  DiffractionPatterns.calc_error_multiangle(probe, rho, data,th_fine_grid(jj),ki_o,kf_o,X,Y,Z);                
            end
            
            figure(fig_num);
            plot(th_fine_grid,err_thin_grid);
            hold on;
            plot(dth_nominal,err_0,'ob');
            cstx = err_0-grad_calc*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
            plot(th_fine_grid,cstx + grad_calc* th_fine_grid);
            title(['theta = ' num2str(dth_nominal)]);
            
            
            
        end
        
        function display_rho_3D(rho,fig_num)
            
            
            % patch( isosurface( smooth3( img, 'gaussian', 7 ) ), 'EdgeColor', 'none', 'FaceColor', 'y' );
            isosurface( ...
                X,Y,Z, ...
                smooth3( abs(NW), 'gaussian', 9 ), 0.5, ...
                smooth3( angle( NW ), 'gaussian', 9 ) );
            axis square; axis tight; axis off; colorbar; lighting gouraud;
            % axis image;
            % axis( [ 0 arr(1) 0 arr(2) 0 arr(3) ] );
            camlight;
            lighting gouraud;
            
        end
    end
end