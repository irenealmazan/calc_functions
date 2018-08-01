classdef SampleGenerator
    % This library contains all the functions to generate the sample
    properties(Constant)
    end
    
    
    methods(Static)
       
        
        function [ imgComplex ] = generateRandomPhase( img, mult,maxphase)
               [ x, y, z ] = meshgrid( 0:size( img, 1 )-1, 0:size( img, 2 )-1, 0:size( img, 3 )-1 );
            x2 = x - mean( x(:) ); x2 = x2 / max( x2(:) );
            y2 = y - mean( y(:) ); y2 = y2 / max( y2(:) );
            z2 = z - mean( z(:) ); z2 = z2 / max( z2(:) );
            temp = rand(3);
            temp = temp'*temp;
            [ R, ~ ] = eig(temp);    % random rotation matrix
            pts = R * [ x2(:) y2(:) z2(:) ]';
            phas = reshape( -pi + 2*pi*sin( mult * 2*pi * ( sum( ( R*pts ).^2 ) ) ), size( img ) );
            phas = maxphase/max( phas(:) ) * phas;
            
            imgComplex = img .* exp( 1i * phas );
            
            if sum(sum(sum(isnan(imgComplex))))
                display('Sample contains NaN')
            end

        end
        
        
        function [ img, fullMap ] = getVoronoiCell( arraySize, pts)
            
            
            
            pts = [ round( arraySize/2 ) ; abs( pts ) ];
            
            %     pts
            [ x, y, z ] = meshgrid( ...
                arraySize(1) * linspace( 0, 1, arraySize(1) ), ...
                arraySize(2) * linspace( 0, 1, arraySize(2) ), ...
                arraySize(3) * linspace( 0, 1, arraySize(3) )  ...
                );
            
            samplePts = [ x(:) y(:) z(:) ];
            
            dist = sqrt( ...
                ( samplePts(:,1) - pts(:,1)' ).^2 + ...
                ( samplePts(:,2) - pts(:,2)' ).^2 + ...
                ( samplePts(:,3) - pts(:,3)' ).^2   ...
                );
            
            [ ~, fullMap ] = min( dist, [], 2 );
            fullMap = reshape( fullMap, arraySize );
            
            img = double( reshape( ...
                ( dist(:,1) == min( dist, [], 2 ) ), ...
                arraySize(1), arraySize(2), arraySize(3) ...
                ) );
            
            [ x, y, z ] = meshgrid( 1:size( img, 1 ), 1:size( img, 2 ), 1:size( img, 3 ) );
            temp = [ x(:) y(:) z(:) ];
            temp = temp( img(:)~=0, : );
            centroid = mean( temp );
            shift = round( -centroid + arraySize/2 );
            img = circshift( img, shift );
            
        end
        
    end
end