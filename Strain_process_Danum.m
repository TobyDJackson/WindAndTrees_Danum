%% The script gives an example of how to filter strain data
% Please do not share without my permission (tobydjackson@gmail.com)


% LOAD UP THE DATA
C1_info=[130.00000	130.00000	130.00000	130.00000	130.00000	130.00000
-10.07000	-12.19000	-11.08000	-14.06000	-10.24000	-10.20000
89.00000	152.00000	138.00000	30.00000	22.00000	270.00000]';
load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\Strain_data_test\Danum\C1_raw.mat')
load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\Strain_data_test\Danum\C1_cup.mat')
load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\Strain_data_test\Danum\C3_sonic.mat')
plot((C1_raw(:,4)))

%% Have a look at the wind data
colors=brewermap(8,'dark2');
peaks_C3sonic=find(C3_sonic(:,2)>3.5); 
histogram(C3_sonic(peaks_C3sonic,3),'FaceColor',colors(3,:),'FaceAlpha',0.5)
[figure_handle,count,speeds,directions,Table] = WindRose(C3_sonic(:,3),(C3_sonic(:,2)),'nspeeds',7,'anglenorth',0,'angleeast',90,'labels',{'N (0°)','S (180°)','E (90°)','W (270°)'},'freqlabelangle',45);

% =======================================
%% 1. Filter and convert from mV data to strain 
% =======================================

%% Part 1. A: Filter and convert from raw data to strain 
data=C1_raw;% % Select a subset of the data to make it quicker.
sample=1:3.3e5;
for col=2:7
    col_data=data(:,col);
    %Filtering for outliers and NaN's
    col_data(find(col_data>=50))=0; %positive outliers
    col_data(find(col_data<=-50))=0; %negative outliers
    col_data(find(isnan(col_data==1)))=0; %NaN's
    %Calibrate with arm length and voltage to extension gradient 'm'
    col_data=col_data/(C1_info(col-1,1)*C1_info(col-1,2));
    if col==2
        data_out=col_data;
    else
        data_out=cat(2,data_out,col_data);% cat together the infor from all the pairs
    end
end % end loop over cols
T_cal=cat(2,data(:,1),data_out);

%% B: Calculate NE, Northward and Eastward facing strain  
for col=2:2:6
    unsmoothed_NE(:,col)=  T_cal(:,col).*cosd(C1_info(col-1,3))+T_cal(:,col+1).*sind(C1_info(col,3));
    unsmoothed_NE(:,col+1)=T_cal(:,col).*sind(C1_info(col-1,3))+T_cal(:,col+1).*cosd(C1_info(col,3));
end
unsmoothed_NE(:,1)=T_cal(:,1);
unsmoothed_NE_sample=unsmoothed_NE(sample,:);
% this version of the data deposited online (EIDC)

% ====================
%% 2. Subtract offsets
% ====================

%% C: Test a running mean vs mode
col=3;
[mode_NE temp_mode] =  Running_mode(unsmoothed_NE(:,col),5000);
mean_NE=unsmoothed_NE(:,col)-running_mean(unsmoothed_NE(:,col),5000);
edges=-5e-4:1e-6:5e-4;
subplot(1,2,1)
histogram(mean_NE(sample),edges);
hold on
xlim([-0.8e-4 0.8e-4])
title('Running mean')
subplot(1,2,2)
histogram(mode_NE(sample),edges);
hold on
xlim([-0.3e-4 0.3e-4])
title('Running mode')

%% C: Smooth using a running mode - test a few different window lengths
% This is done with a smaller sample of data to save time
for col=2:7
    [smoothed1000_NE(:,col-1)   modes1000(:,col-1)]=  Running_mode(unsmoothed_NE_sample(:,col),1000);
    [smoothed5000_NE(:,col-1)   modes5000(:,col-1)]=  Running_mode(unsmoothed_NE_sample(:,col),5000);
    [smoothed15000_NE(:,col-1) modes15000(:,col-1)]=Running_mode(unsmoothed_NE_sample(:,col),15000);
    [smoothed45000_NE(:,col-1) modes45000(:,col-1)]=Running_mode(unsmoothed_NE_sample(:,col),45000);
end

%% D: Calculate MaxStrain, find maxima and estimate CWS
data=cat(2,unsmoothed_NE_sample(:,1),smoothed1000_NE,smoothed5000_NE,smoothed15000_NE,smoothed45000_NE);
T_MaxStrain(:,1)=data(:,1);  
c=1;
for pair=2:2:24
    c=c+1;
    T_MaxStrain(:,c)=sqrt(data(:,pair).^2+data(:,pair+1).^2); 
end

%%  Plotting:
%    A: Plot out the different running modes to se how much variation they are actually subtracting
for col=1:6;
    subplot(4,1,1)
    plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),unsmoothed_NE_sample(:,col+1));
    title('Unsmoothed'); set(gca,'xticklabel',{[]}) ;
    %ylim([-15e-4 15e-4])
    subplot(4,1,2)
    plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),smoothed1000_NE(:,col));
    title('1000 point (4 minute) running mode'); set(gca,'xticklabel',{[]}) ;
    %ylim([-5e-4 5e-4])
    subplot(4,1,3)
    plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),smoothed1000_NE(:,col));
    title('5000 point (21 minute) running mode'); set(gca,'xticklabel',{[]}) ;
    %ylim([-5e-4 5e-4])
    subplot(4,1,4)
    plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),smoothed1000_NE(:,col));
    title('15000 point (62 minute) running mode'); 
    %ylim([-5e-4 5e-4])
    pause
end

%% B: Plot out the different running modes to se how much variation they are actually subtracting
col=4;
plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),unsmoothed_NE_sample(:,col+1));
hold on
plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),modes1000(:,col));
plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),modes5000(:,col));
plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),modes15000(:,col));
plot(datetime(unsmoothed_NE_sample(:,1), 'ConvertFrom', 'datenum'),modes45000(:,col));
h=legend('Strain data','4 minute mode','21 minute mode','62 minute mode','187 minute mode','location','SouthEast');
legend boxoff
title('Running mode window lengths')
ylabel('strain')
set([gca h], 'FontName', 'Helvetica','FontSize', 12)


%% D: Plot two channels and max strain
subplot(3,1,1)
plot(datetime(data(:,1), 'ConvertFrom', 'datenum'),data(:,8));
title('Strain on North facing side of trunk')
ylabel('strain')

subplot(3,1,2)
plot(datetime(data(:,1), 'ConvertFrom', 'datenum'),data(:,9));
title('Strain on East facing side of trunk')
ylabel('strain')

subplot(3,1,3)
plot(datetime(T_MaxStrain(:,1), 'ConvertFrom', 'datenum'),T_MaxStrain(:,5));
title('Max Strain')
ylabel('strain')


%==================
%% 3. Select maxima
%==================

%% Choose one smoothing window for the following sections
clearvars T_MaxStrain data
data=unsmoothed_NE(:,1);  
T_MaxStrain(:,1)=unsmoothed_NE(:,1);  
for col=2:7
    data(:,col)=Running_mode(unsmoothed_NE(:,col),5000);
end
c=1;
for pair=2:2:6
    c=c+1;
    T_MaxStrain(:,c)=sqrt(data(:,pair).^2+data(:,pair+1).^2); 
end

%% Part 3A: Find maxima at specified time intervals - WARNING - This is slow!
minute=6.944445194676518e-04; tenmin=0.006944444496185; hour=0.041666666744277;
wind_in=C1_cup;
[ C1_cup_1min ] = find_variable_maxes( T_MaxStrain, wind_in, minute, 5,220 );
[ C1_cup_10min ] = find_variable_maxes( T_MaxStrain, wind_in, tenmin, 55,2200 );
[ C1_cup_hour ] = find_variable_maxes( T_MaxStrain, wind_in, hour, 330,14000);
scatter(C1_cup_hour(:,2),C1_cup_hour(:,3))

%% Plot the different resolutin maxima
colors=brewermap(8,'Accent');
for col=4:6
    scatter(C1_cup_1min(:,3),C1_cup_1min(:,col),20,colors(2,:),'filled')
    hold on
    scatter(C1_cup_10min(:,3),C1_cup_10min(:,col),20,colors(1,:),'filled')
    scatter(C1_cup_hour(:,3),C1_cup_hour(:,col),20,colors(5,:),'filled')
    legend('1 min', '10 min', 'hour')
    legend boxoff
    ylim([0 0.9e-3])
    title('C1 Wind - Strain relations')
    xlabel('Wind speed (m/s)')
    ylabel('Max strain')
  
    pause
    close all
end


%===========================================
%% 4. Extrapolation up to Critical Wind Speed (CWS)
%===========================================


