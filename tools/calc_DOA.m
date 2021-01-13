function [DOAvec, DOAang] = calc_DOA(loc, arrayCenterPos)
% [DOAvec, DOAang] = calc_DOA(loc, arrayCenterPos)
% calculates DOA vector and DOA angles for given source locations.
%
% IN:
% loc               cartesian coordinates of source locations
% arrayCenterPos    cartesian coordinates of microphone array center
%
% OUT:
% DOAvec            DOA vector
% DOAang            DOA angles (polar, azimuth)


J = size(loc,1);
DOAvec = loc - repmat(arrayCenterPos, J, 1);

ang_pol = zeros(J, 1);
ang_az = zeros(J, 1);

for i = 1:size(DOAvec,1)
    ang_pol = acos(DOAvec(i,3)/norm(DOAvec(i,:)));
    ang_az = atan2(DOAvec(i,2),DOAvec(i,1));
end
    
DOAang = rad2deg([ang_pol ang_az]);

end