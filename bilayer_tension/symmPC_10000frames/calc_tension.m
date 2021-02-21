clear all;

pp=load('press_profile_lateral');
ppmean=get_press(mean(pp));
bins=length(pp(1,:));
slabs=load('press_slabs.txt');
ST=mean(slabs);
plot(-ST*(bins-1)/2:ST:ST*(bins-1)/2,ppmean)

dim =-ST*(bins-1)/2:ST:ST*(bins-1)/2;
bom=1:(bins-1)/2+1;
top=(bins-1)/2+1:bins;
length(top)

tension_top_L=-trapz(1:length(top),ppmean(top)) %bar*A =0.01 mN/m
tension_bot_L=-trapz(1:length(bom),ppmean(bom)) % bar*A
tension_Normal=trapz(1:length(dim),ppmean) % bar*A





