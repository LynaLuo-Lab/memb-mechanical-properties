%-------------------------------------------------------------------------
% LTF Method for leaflet Ka moduli
%
% Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502
%
% Milka Doktorova, March 2019
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Ka_.m
%
% Identifies the relevant thickness for Ka analysis and calculates 
% the leaflet Ka.
%-------------------------------------------------------------------------

function [Ka] = Ka_(sn,a0,temper,flag,snlabel,headCarbon)

warning('off','curvefit:fit:noStartPoint')

cutoff = 10;
binN = 30;

kb = 1.38064853;

fka = fittype('a*x^2+b');
f = fittype('a*x+b');

ncarb = length(sn(1,:));
chain_corr = abs(corr(sn));

x = median(sn(:,1)*ones(1,ncarb)-sn);
y = chain_corr(1,:);

slopes = [];
for i=2:length(x)
    fit1 = fit((x(1:i))',y(1:i)',f);
    slopes = [slopes fit1.a];
end

p = cumsum(slopes(2:end))./[1:length(x)-2];
data = 100*abs((slopes(2:end)-p)./p);

ptdiff = find(data>cutoff);
if length(ptdiff)<2
    Ka = 0;
    disp('Error: Could not identify two regimes in the correlation vs distance data.')
    return
end
pt = ptdiff(1) + 2;

hh = figure;
set(hh,'visible','off');
head = sn(:,pt);
th = (head-min(sn,[],2));

kres = 0.001;
[den,xi,bw]=ksdensity(th,min(th):kres:max(th));
mth = xi(den==max(den));
x = (mth-xi)./mth;
y = den;

threldif = (mth-th)./mth;

rrdata = [];
for rr = 0.07:-0.01:0.05
    h = figure;
    set(h,'visible','off');
    for i=1:10
        realdata = threldif(threldif>=-rr & threldif<=rr);
        realdata = realdata(:);
        if length(realdata)>=100000
            datavec = realdata;
        else
            datavec = repmat(realdata,ceil(100000/length(realdata)),1);
        end
        block = normrnd(datavec,(mth-bw)/mth);    
        q = qqplot(block);
        qx = q.XData;
        qy = q.YData;
        [lf lg lo] = fit(qx',qy',f);
        residuals = lo.residuals;
        qq = qy(abs(residuals)<0.01);
        fq = (max(qq)-min(qq))/(max(qy)-min(qy));
        rrdata = [rrdata; rr fq];
    end
    h.delete
end
rrvec = mean(reshape(rrdata(:,2),[10,3]));
maxind = find(rrvec == max(rrvec));
rr = 0.07 - (maxind-1)*0.01;

indices = find(x>=-rr & x<=rr);
pmf = -2*kb*temper*log(y)/(a0);
fit1 = fit(x(indices)',pmf(indices)',fka);
Ka = fit1.a;

if flag == 1    
    fileID = fopen(strcat('Ka_avgsn_',snlabel,'.txt'),'w');
    fprintf(fileID,'%d mN/m, carbon %d, thickness %4.2f A, fit within range of %4.2f \n',[round(Ka,0) pt+headCarbon-1 round(mth,2) rr]);
    fclose('all');
    
    h = figure();
    subplot(2,2,1)
    hold on
    histogram(th,'BinWidth',0.1,'Normalization','pdf');
    plot(xi,den,'LineWidth',2);
    xlabel('t');
    ylabel('probability');
    title(strcat('C',num2str(pt+1)));
    set(gca,'fontsize',12);
    
    subplot(2,2,2)
    hold on
    x = (mth-xi)./xi;
    y = den;
    pmf = -2*kb*temper*log(y)/(a0);
    reg = -rr:0.001:rr;
    plot(x,pmf);
    plot(reg./(reg+1),fit1.a*(reg./(reg+1)).^2+fit1.b,'LineWidth',2);
    xlim([-0.25 0.25]);
    
    legend('data','fit: ax^2+b','Location','best');
    xlabel('(t0-t)/t');
    ylabel('PMF');
    title(strcat('C',num2str(pt+1)));
    box on
    set(gca,'fontsize',12);

    subplot(2,2,3)
    hold on   
    x = median(sn(:,1)*ones(1,ncarb)-sn);
    y = chain_corr(1,:);
    plot(3:length(x),p,'x','LineWidth',2);
    plot(2:length(x),slopes,'o','LineWidth',2);
    xlim([1 length(x)+1]);
    xlabel('C (from head)');
    ylabel('slope');
    legend('average','current','Location','best');
    set(gca,'fontsize',12);
    box on
    
    subplot(2,2,4)
    hold on
    plot(1:length(x)+1,cutoff*ones(length(x)+1),'--','LineWidth',2);
    plot(3:length(x),100*abs((slopes(2:end)-p)./p),'s','LineWidth',2);
    xlim([1 length(x)+1]);
    xlabel('C (from head)');
    ylabel('slope difference [%]')
    set(gca,'fontsize',12);
    box on
    
    savefig(h,strcat('avgsn_analysis_',snlabel));  
    print(strcat('avgsn_analysis_',snlabel),'-dpng');
    h.delete
    
    % plot distance vs correlation data and fit
    h = figure();
    hold on
    plot(x,y,'o');
    fit1 = fit((x(1:pt-1))',y(1:pt-1)',f);
    yy1 = fit1.a*x+fit1.b;
    plot(x,yy1);
    plot([x(pt) x(pt)],[min(y) max(y)]);
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    grid on
    box on
    xlabel(strcat('distance from head C [A]'));
    ylabel('z correlation');
    title(strcat('K_A=', num2str(round(Ka,0)),' mN/m'));
    set(gca,'fontsize',16);
    savefig(h,strcat('Ka_avgsn_',snlabel));
    print(strcat('Ka_avgsn_',snlabel),'-dpng');
    h.delete
end

return