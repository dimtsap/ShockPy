clc;clf;clear;close all;set(0,'DefaultFigureWindowStyle','docked');warning off
% % Quartz Hugoniot Sandia 2013 Exp Fit
% Function Knudson_exp(wq,x) //a+b up-c exp(-d up)
% 	wave wq
% 	Variable x
% 	return wq[0]+wq[1]*x-wq[2]*x*exp(-wq[3]*x)
% end
wq=[6.278,1.193,2.505,0.3701];
initial_density_Quartz=2.648;
particle_velocity_Quartz=4:0.1:25;
shock_velocity_Quartz=wq(1)+wq(2)*particle_velocity_Quartz-wq(3)*particle_velocity_Quartz.*exp(-wq(4)*particle_velocity_Quartz);
pressure_Quartz=initial_density_Quartz.*shock_velocity_Quartz.*particle_velocity_Quartz;

% % MgO Hugoniot Exp Fit
% Function
% 	wave wm
% 	Variable x
% 	return wm[0]+wm[1]*x+wq[2]*x^2
% end

%% define the states on the Hugoniot where the release path emanate from
wm=[6.6161,1.4111,-0.016277];
initial_density_MgO=3.584;
particle_velocity_MgO=3.5:0.05:10;
shock_velocity_MgO=wm(1)+wm(2)*particle_velocity_MgO+wm(3)*particle_velocity_MgO.^2;
pressure_MgO=initial_density_MgO.*shock_velocity_MgO.*particle_velocity_MgO;
initial_volume_MgO=1/initial_density_MgO;
volume_MgO=initial_volume_MgO .*(shock_velocity_MgO-particle_velocity_MgO) ./shock_velocity_MgO;



%% Reference Hugoniot
upH=linspace(1,10,100);
usH=wm(1)+wm(2)*upH+wm(3)*upH.^2;
vH=initial_volume_MgO.*(usH-upH)./usH;
pH=initial_density_MgO.*usH.*upH;
%% release MgO
gamma=[1.1]
for k=1:length(gamma)
    
    for ii=1:length(pressure_MgO)
        P1=pressure_MgO(ii);
        V1=volume_MgO(ii);
        up1=particle_velocity_MgO(ii);
        us1=shock_velocity_MgO(ii);
        rho0_rho=initial_volume_MgO/V1;
        
        %define volume as the variable
        V=linspace(V1,0.01,1000);
        PH=interp1(vH,pH,V,'cubic');  %interpolation to get the reference Hugoniot Pressure at V
        
        %calculate Es-E0
        
        for jj=1:length(V)
            fun=@(x) (x/V1).^gamma(k) *PH(jj) .*(1-gamma(k)/2.*(initial_volume_MgO./x-1));
            inte(jj)=integral(fun,V1,V(jj));
        end
        
        Es_E0=P1*initial_volume_MgO/2 *(rho0_rho-1)/rho0_rho * (V1./V).^gamma(k) ...
            -(V1./V).^gamma(k) .* inte;
        
        
        Data.Ps(:,ii)=PH .*(1-gamma(k)/2.*(initial_volume_MgO./V-1)) +gamma(k)./V.*(Es_E0);
        Data.ups(:,ii)=up1;
        Vprime=V';
        Data.ups(2:length(V),ii)=real(up1+cumtrapz(V(2:length(V)),sqrt(-diff(Data.Ps(:,ii))./diff(Vprime))));
        
        upIM=linspace(3,17,1000);
        isenP=interp1(Data.ups(2:length(V),ii),Data.Ps(2:length(V),ii),upIM);
        QuartzP=interp1(particle_velocity_Quartz,pressure_Quartz,upIM,'cubic');
        difP=QuartzP-isenP;
        IM.up(ii)=interp1(difP(isnan(difP)==0),upIM(isnan(difP)==0),0);
        IM.UsQ(ii)=interp1(particle_velocity_Quartz,shock_velocity_Quartz,IM.up(ii),'cubic');
        IM.Pq(ii)=interp1(particle_velocity_Quartz,pressure_Quartz,IM.up(ii),'cubic');
        IM.UsMgO(ii)=us1;
        IM.PMgO(ii)=P1;
        
        figure(1);
        hold on;
        plot(V/initial_volume_MgO,Data.Ps(:,ii),'r');
        plot(volume_MgO(ii)/initial_volume_MgO,pressure_MgO(ii),'bo');
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
    plot(particle_velocity_MgO,pressure_MgO,'black');
    plot(particle_velocity_Quartz,pressure_Quartz,'b');
    xlim([3.9,13]);
    ylim([0,700]);
    
    figure(3);
    hold all;
    plot(IM.UsQ,IM.UsMgO,'+');
    xlabel("Us Quartz (km/s)");
    ylabel("Us MgO (km/s)");
    
    header=sprintf('P_MgO  \t Us_MgO \t up_MgO \t  Us_Q  \t P_Q\n ');% GPa \t km/s \t km/s \t Percent  \n');
    IM_Out=[IM.PMgO' IM.UsMgO' IM.up' IM.UsQ' IM.Pq'];
    
%     dlmwrite(['/Users/zixuanye/Documents/MATLAB','MgO_Qtz_IM_Out_',k,'.txt'],header,'delimiter','');
%     dlmwrite(['/Users/zixuanye/Documents/MATLAB','MgO_Qtz_IM_Out_',k,'.txt'],IM_Out,'delimiter','\t','precision','%E','-append');
    
end

%%
Us_mgo_xrd=interp1(IM.UsQ,IM.UsMgO,14.7,'linear');
P_mgo_xrd=interp1(shock_velocity_MgO,pressure_MgO,Us_mgo_xrd,'linear');
up_xrd=interp1(shock_velocity_MgO,particle_velocity_MgO,Us_mgo_xrd);