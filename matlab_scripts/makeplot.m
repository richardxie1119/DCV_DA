%pk_perc_mean = [mean(pks_lr_07(:,1)),mean(pks_lr_05(:,1)),mean(pks_lr_03(:,1)),mean(pks_lr_01(:,1)),mean(pks_lr_005(:,1))];
%pk_err = [std(pks_lr_07(:,1)),std(pks_lr_05(:,1)),std(pks_lr_03(:,1)),std(pks_lr_01(:,1)),std(pks_lr_005(:,1))];
%pk_perc_trunc_mean = [mean(pks_lr_07(:,3)),mean(pks_lr_05(:,3)),mean(pks_lr_03(:,3)),mean(pks_lr_01(:,3)),mean(pks_lr_005(:,3))];
%pk_trunc_err = [std(pks_lr_07(:,3)),std(pks_lr_05(:,3)),std(pks_lr_03(:,3)),std(pks_lr_01(:,3)),std(pks_lr_005(:,3))];
pk_perc_mean2 = [mean(pks_lr_07_2(:,1)),mean(pks_lr_05_2(:,1)),mean(pks_lr_03_2(:,1)),mean(pks_lr_01_2(:,1)),mean(pks_lr_005_2(:,1))];
pk_err2 = [std(pks_lr_07_2(:,1)),std(pks_lr_05_2(:,1)),std(pks_lr_03_2(:,1)),std(pks_lr_01_2(:,1)),std(pks_lr_005_2(:,1))];
pk_perc_trunc_mean2 = [mean(pks_lr_07_2(:,3)),mean(pks_lr_05_2(:,3)),mean(pks_lr_03_2(:,3)),mean(pks_lr_01_2(:,3)),mean(pks_lr_005_2(:,3))];
pk_trunc_err2 = [std(pks_lr_07_2(:,3)),std(pks_lr_05_2(:,3)),std(pks_lr_03_2(:,3)),std(pks_lr_01_2(:,3)),std(pks_lr_005_2(:,3))];

figure();
%errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_mean,pk_err,'LineWidth', 1.5)
hold on
errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_mean2,pk_err2,'LineWidth', 1.5)
%errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_trunc_mean,pk_trunc_err,'LineWidth', 1.5)
errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_trunc_mean2,pk_trunc_err2,'LineWidth', 1.5)
hold off
xlim([0,0.7])
ylim([0,1])
legend('reconstructed','orignial','NumColumns',1)
legend('boxoff')
ylabel('percent of true peaks detected')
xlabel('percent of data used')

figure();
histogram(pks_lr_03_2(:,1));
ylabel('percent of true peaks detected')
xlabel('number of spectrum')
%%
figure();
%pk_corr = [mean(cell2mat(pks_ints_07(:,3))),mean(cell2mat(pks_ints_05(:,3))),mean(cell2mat(pks_ints_03(:,3))),mean(cell2mat(pks_ints_01(:,3))),mean(cell2mat(pks_ints_005(:,3)))];
%pk_corr_std = [std(cell2mat(pks_ints_07(:,3))),std(cell2mat(pks_ints_05(:,3))),std(cell2mat(pks_ints_03(:,3))),std(cell2mat(pks_ints_01(:,3))),std(cell2mat(pks_ints_005(:,3)))];
pk_corr2 = [mean(cell2mat(pks_ints_07_2(:,3))),mean(cell2mat(pks_ints_05_2(:,3))),mean(cell2mat(pks_ints_03_2(:,3))),mean(cell2mat(pks_ints_01_2(:,3))),mean(cell2mat(pks_ints_005_2(:,3)))];
pk_corr_std2 = [std(cell2mat(pks_ints_07_2(:,3))),std(cell2mat(pks_ints_05_2(:,3))),std(cell2mat(pks_ints_03_2(:,3))),std(cell2mat(pks_ints_01_2(:,3))),std(cell2mat(pks_ints_005_2(:,3)))];

%errorbar([0.7,0.5,0.3,0.1,0.05],pk_corr,pk_corr_std,'LineWidth', 1.5)
hold on 
errorbar([0.7,0.5,0.3,0.1,0.05],pk_corr2,pk_corr_std2,'LineWidth', 1.5)

xlim([0,1])
ylim([.98,1])
%legend('dataset 1','dataset 2')
ylabel('mean peak intensity correlation')
xlabel('percent of data used')

figure();
mdl=fitlm(log10(pks_ints_03_2{48,1}(:,2)),log10(pks_ints_03_2{48,1}(:,3)));
plot(mdl)
title('')
xlabel('true peak intensity (log10)')
ylabel('reconstructed peak intensity (log10)')
title(strcat('R-squared = ',num2str(mdl.Rsquared.Adjusted)))

figure();
histogram(log10(cell2mat(pks_ints_03_2(:,4))),20)
xlabel('absolute peak intensity error (log10)')
ylabel('number of peaks')
%% SNR
snr_orig = pks_ints_03_2{48,1}(:,2)/noise_03_2(48,1);
snr_recon = pks_ints_03_2{48,1}(:,3)/noise_03_2(48,2);
figure();
mdl=fitlm(log10(snr_orig),log10(snr_recon));
plot(mdl)
title('')
xlabel('S/N of true peaks')
ylabel('S/N of reconstructed peaks')
title(strcat('R-squared = ',num2str(mdl.Rsquared.Adjusted)))
%%
%err_mean = [mean(ERR_07(:,1)),mean(ERR_05(:,1)),mean(ERR_03(:,1)),mean(ERR_01(:,1)),mean(ERR_005(:,1))];
%err_err = [std(ERR_07(:,1)),std(ERR_05(:,1)),std(ERR_03(:,1)),std(ERR_01(:,1)),std(ERR_005(:,1))];
err_mean2 = [mean(ERR_07_2(:,1)),mean(ERR_05_2(:,1)),mean(ERR_03_2(:,1)),mean(ERR_01_2(:,1)),mean(ERR_005_2(:,1))];
err_err2 = [std(ERR_07_2(:,1)),std(ERR_05_2(:,1)),std(ERR_03_2(:,1)),std(ERR_01_2(:,1)),std(ERR_005_2(:,1))];
figure();
%errorbar([0.7,0.5,0.3,0.1,0.05],err_mean,err_err,'LineWidth', 1.5)
%hold on
errorbar([0.7,0.5,0.3,0.1,0.05],err_mean2,err_err2,'LineWidth', 1.5)
%legend('dataset 1','dataset 2')
ylabel('mean reconstruction error')
xlabel('percent of data used')
%%
pk_perc_mean2 = [mean(pks_lr_07_2(:,2)),mean(pks_lr_05_2(:,2)),mean(pks_lr_03_2(:,2)),mean(pks_lr_01_2(:,2)),mean(pks_lr_005_2(:,2))];
pk_err2 = [std(pks_lr_07_2(:,2)),std(pks_lr_05_2(:,2)),std(pks_lr_03_2(:,2)),std(pks_lr_01_2(:,2)),std(pks_lr_005_2(:,2))];
figure();
%errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_mean,pk_err,'LineWidth', 1.5)
hold on
errorbar([0.7,0.5,0.3,0.1,0.05],pk_perc_mean2,pk_err2,'LineWidth', 1.5)
ylabel('percent of spurious peaks detected')
xlabel('percent of data used')

figure();
histogram(pks_lr_03_2(:,2))
ylabel('number of spectrum')
xlabel('percent of spurious peaks detected')
%%
pks_err = [];
for i=1:length(pks_ints)
    err = abs(pks_ints{i}{1}(:,2)-pks_ints{i}{1}(:,3))./pks_ints{i}{1}(:,2);
    pks_err = [pks_err;err];
end
%%
figure();
histogram(log10(pks_err),20)
xlabel('absolute peak intensity error (%,log10)')
ylabel('number of peaks')
%%
spec = load('spec_out/lung_spec_19.mat');
spec_recon = load('spec_out/lung_spec_recon_19.mat');

figure();
plot(m((200<m)&(m<1000)),spec.varname((200<m)&(m<1000)))
hold on
plot(m((200<m)&(m<1000)),spec_recon.varname((200<m)&(m<1000)))
a=mspeaks(flip(m((500<m)&(m<1000))),flip(spec.varname((500<m)&(m<1000))),'SHOWPLOT',true,'DENOISING',true);
b=mspeaks(flip(m((500<m)&(m<1000))),flip(spec_recon.varname((500<m)&(m<1000))),'SHOWPLOT',true,'DENOISING',true);
%%
ms_file = {'mat_out/pkslr_lung_800_0005.mat',
    'mat_out/pkslr_lung_800_005.mat',
    'mat_out/pkslr_lung_800_01.mat',
    'mat_out/pkslr_lung_800_05.mat'
    };
pk_acc = {};
for i=1:length(ms_file)
    load(ms_file{i})
    pk_acc{i} = pks_lr(:,1);
end
figure();
x = [pk_acc{1};pk_acc{2};pk_acc{3};pk_acc{4}]';
g = [repelem(0.005,length(pk_acc{1})), repelem(0.05,length(pk_acc{2})), repelem(0.1,length(pk_acc{3})),repelem(0.3,length(pk_acc{4}))];
boxplot(x, g,'symbol','')
xlabel('percent of transient used for reconstruction')
ylabel('percent of peaks features recovered per spectrum')
