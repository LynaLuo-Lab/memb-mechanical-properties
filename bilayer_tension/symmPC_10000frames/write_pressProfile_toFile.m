press_tensor_orig = load('run_press');
press_tensor_ew = load('run_ewald');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end) + press_tensor_ew(:,2:end);

press_profile = [];
press_profile_lat = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:)) - pzz(i,:)];
    press_profile_lat = [press_profile_lat; 0.5*(pxx(i,:)+pyy(i,:))];
end

dlmwrite('press_profile_lateral', press_profile_lat,' ');
dlmwrite('press_profile', press_profile,' ');

press_tensor_orig = load('run_press_kin');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end);

press_profile = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:))];% - pzz(i,:)];
end

dlmwrite('press_profile_kin_lateral', press_profile,' ');


press_tensor_orig = load('run_press');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end);

press_profile = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:))];% - pzz(i,:)];
end

dlmwrite('press_profile_noew_lateral', press_profile,' ');


press_tensor_orig = load('run_press_bond');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end);

press_profile = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:))];% - pzz(i,:)];
end

dlmwrite('press_profile_bond_lateral', press_profile,' ');


press_tensor_orig = load('run_press_nonbond');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end);

press_profile = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:))];% - pzz(i,:)];
end

dlmwrite('press_profile_nonbond_lateral', press_profile,' ');

% figure
press_tensor_orig = load('run_ewald');

num_slices = length(press_tensor_orig(1,2:end))/3;
num_frames = length(press_tensor_orig(:,1));

press_tensor = press_tensor_orig(:,2:end);

press_profile = [];
pxx = press_tensor(:,1:3:end);
pyy = press_tensor(:,2:3:end);
pzz = press_tensor(:,3:3:end);

for i = 1:num_frames
    press_profile = [press_profile; 0.5*(pxx(i,:)+pyy(i,:))];% - pzz(i,:)];
end

dlmwrite('press_profile_ewald_lateral', press_profile,' ');
