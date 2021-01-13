function [loc, Delta_t_i] = gen_searchGrid(micPos, dim1, dim2, mode, c)
% [loc, Delta_t_i] = gen_searchGrid(micPos, dim1, dim2, mode, c)
% converts grid points to location coordinates (near field)/DOA vectors(far
% field) and computes TDOAs
%
% IN:
% micPos         microphone positions - channels x coordinates
% dim1           grid points first dimension (cartesian mode = x-coordinate; polar mode = distance; spherical mode = polar angle)
% dim2           grid points second dimension (cartesian mode = y-coordinate; polar mode = angle; spherical mode = azimuth angle)
% c              speed of sound
%
% OUT:
% loc            location coordinates or DOA vectors -  candidate locations x microphone pairs
% Delta_t_i      TDOAs - candidate locations x microphone pairs


% number of microphones
M = size(micPos, 1);
% nubmer of microphone pairs
P = M*(M-1)/2;

N_dim1 = length(dim1);
N_dim2 = length(dim2);

    
i = 0;
for n_dim1 = 1:N_dim1
    for n_dim2 = 1:N_dim2
        
        i = i + 1;
        switch mode
            case 'cartesian'
                x = dim1(n_dim1);
                y = dim2(n_dim2);
                % location in cartesian coordinates
                loc(i, :) = [x y];
                % TDOA
                Delta_t_i(i,:) = loc2TDOA(micPos, [x y], c);
            case 'polar'
                r = dim1(n_dim1);
                theta = dim2(n_dim2);
                x = r*cos(deg2rad(theta));
                y = r*sin(deg2rad(theta));
                % location in cartesian coordinates
                loc(i, :) = [x y];
                % TDOA
                Delta_t_i(i,:) = loc2TDOA(micPos, [x y], c);
            case 'spherical'
                ang_pol = dim1(n_dim1);
                ang_az = dim2(n_dim2);
                % location in far field spherical coordinates
                [Delta_t_i(i,:), loc(i,:)] = spher2TDOA(micPos, ang_pol,  ang_az, c);
                if ang_pol == 0 || ang_pol == 180
                    break % we do not need to loop further over ang_az because DOAs won't change
                end
        end
    end
end
    
    

end


function delta_t = loc2TDOA(micPos, loc, c)


M = size(micPos,1);
P = M*(M-1)/2;

% propagation time
dist = sqrt(sum((micPos - repmat(loc, size(micPos, 1), 1)).^2, 2));
propTime = dist/c;

delta_t = zeros(P,1);
p = 0;
for mprime = 1:M
    for m = mprime+1:M
        p = p+1;
        delta_t(p) = propTime(m) - propTime(mprime);
    end
end

end


function [delta_t, DOA_vec] = spher2TDOA(micPos, ang_pol, ang_az, c)


M = size(micPos,1);
P = M*(M-1)/2;

DOA_vec = [sin(deg2rad(ang_pol))*cos(deg2rad(ang_az)),...
    sin(deg2rad(ang_pol))*sin(deg2rad(ang_az)),...
    cos(deg2rad(ang_pol))];

b = -transp(DOA_vec);

delta_t = zeros(P,1);
p = 0;
for mprime = 1:M
    for m = mprime+1:M
        p = p+1;
        
        a = micPos(m,:).' - micPos(mprime,:).';

        delta_t(p) = (norm(a)/c)*(a.'*b/(norm(a)*norm(b)));

    end
end

end
