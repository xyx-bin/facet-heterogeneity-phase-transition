%function [y1,...,yN] = myfun(x1,...,xM) .��Ϊ myfun �ĺ�����
%�ú����������� x1,...,xM ��������� y1,...,yN
function dx=myfun(t,x)
%�޸ı���Ϊȫ�ֱ���
global beta1 beta2 k1 k2 gamma
dx=zeros(1,1);
dx(1)=-gamma*x+beta1*k1*x*(1-x)+beta2*k2*x^2*(1-x);
end
