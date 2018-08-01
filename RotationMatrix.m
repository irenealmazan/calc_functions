classdef RotationMatrix
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
       
         function [Ry_rock, Ry_rock_deriv] = rock_curve(dth,dth_perturb)
            
            % This function calculates the rotation matrix which turns the qbragg peak in order to scan the
            % rocking curve and its derivative with respect to dth. Needs to be adapted for each experimental
            % set-up
            
            Ry_rock = [cosd(dth) 0 sind(dth);
                0 1 0;
                -sind(dth) 0 cosd(dth)];
            
%             Ry_rock_perturb = [cosd(dth_perturb) 0 sind(dth_perturb);
%                 0 1 0;
%                 -sind(dth_perturb) 0 cosd(dth_perturb)];
%             
            Ry_rock_perturb_deriv = [-sind(dth_perturb) 0 +cosd(dth_perturb);
                0 0 0;
                -cosd(dth_perturb) 0 -sind(dth_perturb)];
            
             Ry_rock_deriv = (pi/180)*Ry_rock*Ry_rock_perturb_deriv;
              
        end
        
        function [Ry,Rx] = detector(del,gam)
            
            % This function writes the rotational matrix which describe the position
            % of the detector with respect to the outgoing beam. The reference frame is
            % sector 26 of APS
            
            Ry = [cosd(del) 0 sind(del);
                0 1 0;
                -sind(del) 0 cosd(del)];
            
            Rx=[1 0 0;
                0 cosd(gam) -sind(gam);
                0 sind(gam) cosd(gam)];
        end
        
        
        function [Rz,Ry_th,Ry_del,Rx] = NW(phi,th,del,gam)
            
            % This function writes the rotational matrix which describe the position
            % of the detector with respect to the outgoing beam. The reference frame is
            % sector 26 of APS
            
            Rz=[cosd(phi) -sind(phi) 0;
                sind(phi) cosd(phi) 0;
                0 0 1];
            
            Ry_th=[cosd(th) 0 sind(th);
                0 1 0;
                -sind(th) 0 cosd(th)];
            
            Ry_del = [cosd(-del) 0 sind(-del);
                0 1 0;
                -sind(-del) 0 cosd(-del)];
            
            Rx=[1 0 0;
                0 cosd(gam) -sind(gam);
                0 sind(gam) cosd(gam)];
        end
        
        
    end
    

end