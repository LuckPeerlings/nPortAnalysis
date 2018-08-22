close all
clear all

figure
AX1 = subplot(2,1,1); hold
AX2 = subplot(2,1,2); hold

figure
AX3 = subplot(2,1,1); hold
AX4 = subplot(2,1,2); hold

[FILE,DIR]=uigetfile;
Z1 = load([DIR,FILE]);


for i=1:4
    C_fitted_1(:,i) = polyval(Z1.PolReal(i,:),Z1.f_cal,[],Z1.MU_Real(i,:)) + 1i.*polyval(Z1.PolImag(i,:),Z1.f_cal,[],Z1.MU_Imag(i,:));
end

[FILE,DIR]=uigetfile;
Z2 = load([DIR,FILE]);

for i=1:4
    C_fitted_2(:,i) = polyval(Z2.PolReal(i,:),Z2.f_cal,[],Z2.MU_Real(i,:)) + 1i.*polyval(Z2.PolImag(i,:),Z2.f_cal,[],Z2.MU_Imag(i,:));
end

% plot(AX1,Z1.f_cal,abs(C_fitted_1))
% plot(AX2,Z1.f_cal,angle(C_fitted_1))
% plot(AX1,Z2.f_cal,abs(C_fitted_2))
% plot(AX2,Z2.f_cal,angle(C_fitted_2))


plot(AX3,Z1.f_cal,abs(Z1.C(:,1:end)))
plot(AX4,Z1.f_cal,angle(Z1.C(:,1:end)))

plot(AX3,Z2.f_cal,abs(Z2.C),'-.')
plot(AX4,Z2.f_cal,angle(Z2.C),'-.')