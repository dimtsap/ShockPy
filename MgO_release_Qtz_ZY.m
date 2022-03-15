clc;clf;clear;close all;set(0,'DefaultFigureWindowStyle','docked');warning off
% % Quartz Hugoniot Sandia 2013 Exp Fit
% Function Knudson_exp(wq,x) //a+b up-c exp(-d up)
% 	wave wq
% 	Variable x
% 	return wq[0]+wq[1]*x-wq[2]*x*exp(-wq[3]*x)
% end
wq=[6.278,1.193,2.505,0.3701];
rho0q=2.648;
upq=4:0.1:25;
usq=wq(1)+wq(2)*upq-wq(3)*upq.*exp(-wq(4)*upq);
pq=rho0q.*usq.*upq;

% % MgO Hugoniot Exp Fit
% Function
% 	wave wm
% 	Variable x
% 	return wm[0]+wm[1]*x+wq[2]*x^2
% end

%% define the states on the Hugoniot where the release path emanate from
wm=[6.6161,1.4111,-0.016277];
rho0m=3.584;
upm=3.5:0.05:10;
usm=wm(1)+wm(2)*upm+wm(3)*upm.^2;
pm=rho0m.*usm.*upm;
v0m=1/rho0m;
vm=v0m .*(usm-upm) ./usm;



%% Reference Hugoniot
upH=linspace(1,10,100);
usH=wm(1)+wm(2)*upH+wm(3)*upH.^2;
vH=v0m.*(usH-upH)./usH;
pH=rho0m.*usH.*upH;
%% release MgO
gamma=[1.1]
for k=1:length(gamma)
    
    for ii=1:length(pm)
        P1=pm(ii);
        V1=vm(ii);
        up1=upm(ii);
        us1=usm(ii);
        ita=v0m/V1;
        
        %define volume as the variable
        V=linspace(V1,0.244,1000);
        PH=interp1(vH,pH,V,'cubic');  %interpolation to get the reference Hugoniot Pressure at V
        
        %calculate Es-E0
        
        for jj=1:length(V)
            fun=@(x) (x/V1).^gamma(k) *PH(jj) .*(1-gamma(k)/2.*(v0m./x-1));
            inte(jj)=integral(fun,V1,V(jj));
        end
        
        Es_E0=P1*v0m/2 *(ita-1)/ita * (V1./V).^gamma(k) ...
            -(V1./V).^gamma(k) .* inte;
        Data.Ps(:,ii)=PH .*(1-gamma(k)/2.*(v0m./V-1)) +gamma(k)./V.*(Es_E0);
        Data.ups(:,ii)=up1;
        Vprime=V';
        Data.ups(2:length(V),ii)=real(up1+cumtrapz(V(2:length(V)),sqrt(-diff(Data.Ps(:,ii))./diff(Vprime))));
        
        upIM=linspace(3,17,1000);
        isenP=interp1(Data.ups(2:length(V),ii),Data.Ps(2:length(V),ii),upIM);
        QuartzP=interp1(upq,pq,upIM,'cubic');
        difP=QuartzP-isenP;
        IM.up(ii)=interp1(difP(isnan(difP)==0),upIM(isnan(difP)==0),0);
        IM.UsQ(ii)=interp1(upq,usq,IM.up(ii),'cubic');
        IM.Pq(ii)=interp1(upq,pq,IM.up(ii),'cubic');
        IM.UsMgO(ii)=us1;
        IM.PMgO(ii)=P1;
        
        figure(1);
        hold on;
        plot(V/v0m,Data.Ps(:,ii),'r');
        plot(vm(ii)/v0m,pm(ii),'bo');
        xlabel("V/V0 ");
        ylabel("Pressure (GPa)");
        
        figure(2);
        hold all;
        plot(Data.ups(:,ii),Data.Ps(:,ii),'red');
        xlabel("up (km/s)");
        ylabel("Pressure (GPa)");
        
    end
    
    figure(2);
    hold all;
    plot(upm,pm,'black');
    plot(upq,pq,'b');
    xlim([3.9,13]);
    ylim([0,700]);
    
    figure(3);
    hold all;
    plot(IM.UsQ,IM.UsMgO,'+');
    xlabel("Us Quartz (km/s)");
    ylabel("Us MgO (km/s)");
    
    header=sprintf('P_MgO  \t Us_MgO \t up_MgO \t  Us_Q  \t P_Q\n ');% GPa \t km/s \t km/s \t Percent  \n');
    IM_Out=[IM.PMgO' IM.UsMgO' IM.up' IM.UsQ' IM.Pq'];
    
    dlmwrite(['/Users/zixuanye/Documents/MATLAB','MgO_Qtz_IM_Out_',k,'.txt'],header,'delimiter','');
    dlmwrite(['/Users/zixuanye/Documents/MATLAB','MgO_Qtz_IM_Out_',k,'.txt'],IM_Out,'delimiter','\t','precision','%E','-append');
    
end

%%
Us_mgo_xrd=interp1(IM.UsQ,IM.UsMgO,14.7,'linear');
P_mgo_xrd=interp1(usm,pm,Us_mgo_xrd,'linear');
up_xrd=interp1(usm,upm,Us_mgo_xrd);