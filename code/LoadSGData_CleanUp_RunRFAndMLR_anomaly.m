%LoadSGData_CleanUp_RunRFAndMLR_anomaly.m
%For the glider file selected at top of script, this code cleans up the
%relevant CTD and O2 files amd grabs a glider down-dive collected every 5
%days (e.g., at the Argo sampling resolution). It then constructs two 
%types of regression models (random forest and multiple linear regression)
%that relate latitude, longitude, pressure, temperature, and salinity to
%predict oxygen using these "Argo profiles". The predictions from the two
%methods are compared with the measured oxygen values at all reaming glider
%dives to assess performance.
%pdlogan@uw.edu Updated: 01/14/21

%% Choose glider files
% fname1 = 'CTD_659.nc'; fname2 = 'O2_659.nc'; fname3 = 'RF_MLR_O2_ResultsAtTestPts_659.nc';
fname1 = 'CTD_660.nc'; fname2 = 'O2_660.nc'; fname3 = 'RF_MLR_O2_ResultsAtTestPts_660.nc';

% anom = 0; %then run regressions using oxygen values
anom = 1; %then run regressions using oxygen ANOMALY values

%% Load CTD data
p = ncread(fname1,'pressure'); %dbar
lat = ncread(fname1,'latitude');
lon = ncread(fname1,'longitude');
t = ncread(fname1,'temperature'); %Conservative Temperature
s = ncread(fname1,'salinity'); %Absolute Salinity
divenum = ncread(fname1,'dives');
time = ncread(fname1,'time'); %seconds since 1970-01-01
secperday = 86400;
timeday = time/secperday;
ctdstarttimeday = timeday(1);
timeday = timeday-ctdstarttimeday;

badp = find(isnan(p)==1);
badlat = find(isnan(lat)==1);
badlon = find(isnan(lon)==1);
badt = find(isnan(t)==1);
bads = find(isnan(s)==1);
badtime = find(isnan(time)==1);
nanidx = unique([badp;badlat;badlon;badt;bads;badtime]);

timeday(nanidx) = []; time(nanidx) = []; p(nanidx) = []; lat(nanidx) = [];
lon(nanidx) = []; t(nanidx) = []; s(nanidx) = []; divenum(nanidx) = [];
clear badp badlat badlon badt bads badtime nanidx

npts = length(p);
divenumvect = unique(divenum);
ndivelegs = length(divenumvect);

%% Plot all temperature data and find start day of each dive
startday = [];
for i = 1:ndivelegs
    idx = find(divenum==divenumvect(i));
%     scatter(time(idx),p(idx),[],t(idx),'filled'); colorbar; hold on;
    startday = [startday;timeday(idx(1))];
end
% set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure')

%grab start of down dives only
startday_down = startday(1:2:end); %see that last dive started at 86.2703 days

%% Load in oxygen file and remove NaNs
p_ox = ncread(fname2,'pressure'); %dbar
ox = ncread(fname2,'oxygen'); %micromoles/L
time_ox = ncread(fname2,'time'); %seconds since 1970-01-01
timeday_ox = time_ox/secperday;
timeday_ox = timeday_ox-ctdstarttimeday; %reference relative to first timeday in CTD file nor ox file

badp = find(isnan(p_ox)==1);
badox = find(isnan(ox)==1);
badtime = find(isnan(time_ox));
nanidx = unique([badp;badox;badtime]);
p_ox(nanidx) = []; ox(nanidx) = [];
time_ox(nanidx) = []; timeday_ox(nanidx) = [];
clear badp badox badtime nanidx

npts_ox = length(p_ox);

%% Make a divenum variable for oxygen
% divestartidx_ox = NaN*ones(length(startday_down),1);
divenum_ox = NaN*ones(npts_ox,1);
for i = 1:length(startday) %loop over each dive
%     divestartidx_ox = find(startday_down(i)==timed_ox);
    if i == 1
        indive_idx = find(timeday_ox<startday(i+1));
        divenum_ox(indive_idx) = divenumvect(i);
%         figure(2); scatter(timeday_ox(indive_idx),p_ox(indive_idx),[],ox(indive_idx),'filled'); colorbar;
%         hold on; scatter(timeday(divenum==divenumvect(i)),p(divenum==divenumvect(i)));
%         keyboard()
    elseif i == length(startday)
        indive_idx = find(timeday_ox>=startday(i));
        divenum_ox(indive_idx) = divenumvect(i);
%         figure(2); scatter(timeday_ox(indive_idx),p_ox(indive_idx),[],ox(indive_idx),'filled'); colorbar;
%         hold on; scatter(timeday(divenum==divenumvect(i)),p(divenum==divenumvect(i)));
%         keyboard()
    else
        indive_idx = find(timeday_ox>=startday(i)&timeday_ox<startday(i+1));
        divenum_ox(indive_idx) = divenumvect(i);
%         figure(2); scatter(timeday_ox(indive_idx),p_ox(indive_idx),[],ox(indive_idx),'filled'); colorbar;
%         hold on; scatter(timeday(divenum==divenumvect(i)),p(divenum==divenumvect(i)));
%         keyboard()
    end
end

%% For each dive leg linearly interpolate oxygen data to same p/time as t/s
diveswithoutox = [];
ox_interp = [];
for i = 1:ndivelegs
    thisdive = find(divenum==divenumvect(i));
    thisdive_ox = find(divenum_ox==divenumvect(i));
    if isempty(thisdive_ox)==0&length(thisdive)>1&length(thisdive_ox)>1
        %need this check b/c some dives lack oxygen data & some dives only
        %have a single point in 659 file
        if length(unique(p_ox(thisdive_ox)))==length(p_ox(thisdive_ox))
            ox_tmp = interp1(p_ox(thisdive_ox),ox(thisdive_ox),p(thisdive)); %interpolate to CTD p-levels
        elseif length(unique(time_ox(thisdive_ox)))==length(time_ox(thisdive_ox)) %only do if repeated pressure levels would make prev interp fail (e.g., near-surf)
            ox_tmp = interp1(time_ox(thisdive_ox),ox(thisdive_ox),time(thisdive)); %interpolate to CTD time-levels
        else
            disp('Debug.'); keyboard()
        end
%         figure(3); subplot(1,2,1); scatter(time_ox(thisdive_ox),p_ox(thisdive_ox),[],ox(thisdive_ox),'filled'); colorbar;
%         set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure'); title('Oxygen Profile (Measured)'); colorbar
%         subplot(1,2,2); scatter(time(thisdive),p(thisdive),[],ox_tmp,'filled'); title('Oxygen Profile (Interpolated to T/S Levels)'); colorbar
%         set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure')
%         title('Oxygen Profile (Interpolated to T/S Levels)');
        ox_interp = [ox_interp;ox_tmp];
    else %then dive leg lacks oxygen data
        diveswithoutox = [diveswithoutox;divenumvect(i)];
    end
end

%% Remove dives without oxygen data
nooxind = [];
for i = 1:length(diveswithoutox)
    nooxind = [nooxind;find(divenum==diveswithoutox(i))];
end
%Make sure grabbed all nooxind
if length(nooxind)~= (length(divenum)-length(ox_interp))
    disp('Debug.'); keyboard()
end

%Remove points without oxygen
p(nooxind) = []; lat(nooxind) = []; lon(nooxind) = [];
t(nooxind) = []; s(nooxind) = []; divenum(nooxind) = [];
time(nooxind) = []; timeday(nooxind) = [];
npts = length(divenum);
ndivelegs = length(unique(divenum));
divenumvect = unique(divenum);
% ox_interp is already the correct size


%% Divide test vs training data
%grab start times of profiles approx every 5 days
tmp = 0:5:startday_down(end);
argoprof_diveidx = NaN*ones(length(tmp),1);
for i = 1:length(tmp)
    argoprof_diveidx(i) = find(abs(startday_down-tmp(i))==min(abs(startday_down-tmp(i))));
end

fprintf('There are %d "Argo profiles" in this glider dataset \n',length(tmp));
%since only 86 days of data, this leads to only 18 Argo "profiles" at
%5-day intervals

apdi = argoprof_diveidx; clear argoprof_diveidx

train_i = [];
for i = 1:length(apdi)
    train_i = [train_i;find(divenum==apdi(i))];
end
%24,797/1,341,072=1.85% data points are training data for SG 659
%20,955/1,273,697=1.65% data points are training data for SG 660
test_i = [1:npts]'; test_i(train_i)=[];


%% Separate out train and test dat then plot
figure; scatter(time(train_i),p(train_i),[],ox_interp(train_i),'filled');
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');
title('"Argo Profiles" of Oxygen [Train]'); colorbar; caxis([175 400]); ylim([-5 1050])
figure; scatter(time(test_i),p(test_i),[],ox_interp(test_i),'filled');
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');
title('Glider Dives of Oxygen [Test]'); colorbar; caxis([175 400]);
ylim([-5 1050]); xlim([1.5565e9 1.5643e9])

%% De-mean oxygen data
if anom == 1
    %Interpolate each oxygen profile to standard levels
    tmp = p_ox;
    tmptmp = diff(tmp);
    figure; scatter(1:length(tmptmp),tmptmp)
    ylim([0 10])
    
    tmp = p(train_i); tmptmp = diff(tmp);
    figure; scatter(1:length(tmptmp),tmptmp)
    fprintf('Median spacing of oxygen measurements = %f dbar \n',median(tmptmp))
    clear tmp tmptmp
    disp('We will interpolate them to 2 dbar for now though')
    disp('Note that currently we are again interpolating the oxygen values which were previously interpolated once already to T/S levels');
    
    interp_int = 2;
    min_p = ceil(min(p(train_i)));
    max_p = floor(max(p(train_i)));
    ndivelegs_train = length(unique(divenum(train_i)));
    p_interpstd = min_p:interp_int:max_p; %standard levels to interpolate to
    ox_interpstd = NaN*ones(length(p_interpstd),ndivelegs_train);
    for i = 1:ndivelegs_train
        thisdive = find(divenum==apdi(i));
        %This interp1 function does not extrapolate so extra p-levls will be
        %filled with NaNs
        ox_tmp = interp1(p(thisdive),ox_interp(thisdive),p_interpstd);
        %     figure(13); subplot(1,2,1); scatter(ox_interp(thisdive),p(thisdive))
        %     set(gca,'YDir','reverse'); xlabel('Oxygen'); ylabel('Pressure');
        %     subplot(1,2,2); scatter(ox_tmp,p_interpstd);
        %     set(gca,'YDir','reverse'); xlabel('Oxygen'); ylabel('Pressure');
        ox_interpstd(:,i) = ox_tmp;
    end
    
    %Find average profile across all oxygen profiles
    ox_avg = nanmean(ox_interpstd,2);
    p_avg = p_interpstd';
    p_avg(find(isnan(ox_avg)==1))=[];
    ox_avg(find(isnan(ox_avg)==1))=[];
    figure; scatter(ox_avg,p_avg)
    set(gca,'YDir','reverse'); xlabel('Oxygen'); ylabel('Pressure');
    
    %Interpolate this average profile to the depth levels of each profile and
    %then remove it from each profile (did for training and test all at once)
    ox_anom = [];
    ox_avg_interp = [];
    for i = 1:ndivelegs
        thisdive = find(divenum==divenumvect(i));
        ox_tmp = interp1(p_avg,ox_avg,p(thisdive));
        ox_tmp2 = ox_interp(thisdive)-ox_tmp;
        ox_anom = [ox_anom;ox_tmp2];
        ox_avg_interp = [ox_avg_interp;ox_tmp];
        %     figure(15); clf; scatter(ox_interp(thisdive),p(thisdive))
        %     set(gca,'YDir','reverse'); xlabel('Oxygen'); ylabel('Pressure');
        %     hold on; scatter(ox_avg,p_avg)
        %     figure(16); clf; scatter(ox_tmp2,p(thisdive))
        %     set(gca,'YDir','reverse'); xlabel('Oxygen Anomaly'); ylabel('Pressure');
        %     keyboard()
    end
end

%% Predict oxygen using random forest regression (with time as a predictor)
%then assess performance
%Does significantly worse than case without time (assessed using median
% absolute error and by looking at the histogram of the absolute errors

%Set-up training/test arrays
% X_train = [lat(train_i),lon(train_i),time(train_i),p(train_i),t(train_i),s(train_i)]; Y_train = [ox_interp(train_i)];
% X_test = [lat(test_i),lon(test_i),time(test_i),p(test_i),t(test_i),s(test_i)]; Y_test = [ox_interp(test_i)];

%Construct random forest model, predict for glider test data, and calculate
%absolute error relative to test observations
% rng(1); % For reproducibility
% Mdl = TreeBagger(100,X_train,Y_train,'Method','regression','OOBPredictorImportance','On');
% Y_pred = predict(Mdl,X_test);
% 
% AE = Y_pred-Y_test;
% 
% figure; scatter(1:length(Y_test),AE);
% figure; histogram(AE);
% figure; scatter(1:length(Y_test),sort(AE));
% fprintf('Median absolute error from RF prediction: \n %f \n',nanmedian(AE));
% 
% figure
% bar(Mdl.OOBPermutedPredictorDeltaError)
% xlabel('Feature Number') 
% ylabel('Out-of-Bag Feature Importance')
% xticklabels({'lat','lon','time','p','t','s'})

%% Predict oxygen using random forest regression (without time as a 
%predictor) then assess performance
%Set-up training/test arrays
if anom == 0
    X_train = [lat(train_i),lon(train_i),p(train_i),t(train_i),s(train_i)]; Y_train = [ox_interp(train_i)];
    X_test = [lat(test_i),lon(test_i),p(test_i),t(test_i),s(test_i)]; Y_test = [ox_interp(test_i)];
elseif anom == 1
    X_train = [lat(train_i),lon(train_i),p(train_i),t(train_i),s(train_i)]; Y_train = [ox_anom(train_i)];
    X_test = [lat(test_i),lon(test_i),p(test_i),t(test_i),s(test_i)]; Y_test = [ox_anom(test_i)];
end

%Construct random forest model, predict for glider test data, and calculate
%absolute error relative to test observations
rng(1); % For reproducibility
nTrees = 100;
Mdl = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPredictorImportance','On');
Y_pred = predict(Mdl,X_test);
AE = Y_pred-Y_test;
fprintf('Median absolute error from RF prediction: \n %f \n',nanmedian(AE));
fprintf('Median absolute absolute error from RF prediction: \n %f \n',nanmedian(abs(AE)));

% figure; scatter(1:length(Y_test),AE);
figure; histogram(AE); xlabel('Absolute Error'); ylabel('Number of Points');
title('RF');
% figure; scatter(1:length(Y_test),sort(AE));

figure
bar(Mdl.OOBPermutedPredictorDeltaError)
xlabel('Feature Number') 
ylabel('Out-of-Bag Feature Importance')
xticklabels({'lat','lon','p','t','s'})

% Figure to explore the appropriate number of trees
% figure
% plot(oobError(Mdl))
% xlabel('Number of Grown Trees')
% ylabel('Out-of-Bag Mean Squared Error')
%could explore different minleafsizes as well later

figure; scatter(time(test_i),p(test_i),[],AE,'filled');
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');
if anom==0
    title('Absolute Error (from RF) for Test Glider Dives of Oxygen'); colorbar;
    clmap(23); caxis([-80 80]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
elseif anom==1
    title('Absolute Error (from RF) for Test Glider Dives of Oxygen Anomaly'); colorbar;
    clmap(23); caxis([-80 80]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
end

%% Predict oxygen using basic multiple linear regression then assess performance
X = [ones(size(Y_train)) X_train(:,1) X_train(:,2) X_train(:,3) X_train(:,4) X_train(:,5)];
[b,~,~,~,stats] = regress(Y_train,X); R2 = stats(1);
Y_pred_mlr = b(1) + b(2)*X_test(:,1) + b(3)*X_test(:,2) + b(4)*X_test(:,3)+ b(5)*X_test(:,4)+ b(6)*X_test(:,5);
AE_mlr = Y_pred_mlr-Y_test;
fprintf('R^2 for MLR fit: \n %f \n',R2);
fprintf('Median absolute error from MLR prediction: \n %f \n',nanmedian(AE_mlr));
fprintf('Median absolute absolute error from MLR prediction: \n %f \n',nanmedian(abs(AE_mlr)));

figure; histogram(AE_mlr); xlabel('Absolute Error'); ylabel('Number of Points');
title('MLR');

figure; scatter(time(test_i),p(test_i),[],AE_mlr,'filled');
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');
if anom==0
    title('Absolute Error (from MLR) for Test Glider Dives of Oxygen'); colorbar;
    clmap(23); caxis([-80 80]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
elseif anom==1
    title('Absolute Error (from MLR) for Test Glider Dives of Oxygen Anomaly'); colorbar;
    clmap(23); caxis([-80 80]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
end

%% Plot measured and predicted oxygen values along "test" glider tracks
figure; scatter(time(test_i),p(test_i),[],Y_test);
if anom==0
    title('Measured Oxygen Values')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
elseif anom==1
    title('Measured Oxygen Anomaly Values'); clmap(23)
    colorbar; caxis([-60 60]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
end
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');

figure; scatter(time(test_i),p(test_i),[],Y_pred);
if anom==0
    title('Predicted Oxygen Values (RF)')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
elseif anom==1
    title('Predicted Oxygen Anomaly Values (RF)'); clmap(23)
    colorbar; caxis([-60 60]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
end
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');

figure; scatter(time(test_i),p(test_i),[],Y_pred_mlr);
if anom==0
    title('Predicted Oxygen Values (MLR)')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
elseif anom==1
    title('Predicted Oxygen Anomaly Values (MLR)'); clmap(23)
    colorbar; caxis([-60 60]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
end
set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');

%% Plot measured and predicted oxygen values along "test" glider tracks
if anom ==1
    figure; scatter(time(test_i),p(test_i),[],Y_test+ox_avg_interp(test_i));
    title('Measured Oxygen Values')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
    set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');

    figure; scatter(time(test_i),p(test_i),[],Y_pred+ox_avg_interp(test_i));
    title('Predicted Oxygen Values (RF)')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
    set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');

    figure; scatter(time(test_i),p(test_i),[],Y_pred_mlr+ox_avg_interp(test_i));
    title('Predicted Oxygen Values (MLR)')
    colorbar; caxis([175 400]); ylim([-5 1050]); xlim([1.5565e9 1.5643e9])
    set(gca,'YDir','reverse'); xlabel('Time'); ylabel('Pressure');
end

%% Write results to netcdf file
% nccreate(fname3,'time','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'pressure','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'ox_measured','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'ox_predicted_RF','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'ox_predicted_MLR','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'absolute_error_RF','Dimensions',{'npts_test',length(test_i),'1',1})
% nccreate(fname3,'absolute_error_MLR','Dimensions',{'npts_test',length(test_i),'1',1})
% 
% ncwrite(fname3,'time',time(test_i))
% ncwrite(fname3,'pressure',p(test_i))
% ncwrite(fname3,'ox_measured',Y_test)
% ncwrite(fname3,'ox_predicted_RF',Y_pred)
% ncwrite(fname3,'ox_predicted_MLR',Y_pred_mlr)
% ncwrite(fname3,'absolute_error_RF',AE)
% ncwrite(fname3,'absolute_error_MLR',AE_mlr)
% 
% ncdisp(fname3)