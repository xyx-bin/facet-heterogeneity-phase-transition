clear;
clc;
syms gamma k1 k2 N average_2;
gamma=0.6;N=500;average_2=2;
yy=[];
Beta1=[0.03:0.01:0.3];%0.3限制k^2/k > 1    0.03限制k^2/k < 10
for i=1:length(Beta1)
    beta1=Beta1(i);
    y = (gamma)/(2*beta1);
    yy(i)=y;
end
plot(yy,Beta1,'linewidth',2);
xlabel('k^2/k');
ylabel('\beta_1');
hold on