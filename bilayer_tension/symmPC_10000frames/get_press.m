function [press] = get_press(profile)

slabs = length(profile);
latT = trapz(1:slabs,profile);
pz=(latT/((slabs-1)))*ones(1,slabs);

press = profile-pz;

end