%function [y1,...,yN] = myfun(x1,...,xM) .名为 myfun 的函数，
%该函数接受输入 x1,...,xM 并返回输出 y1,...,yN
function dx=myfun(t,x)
%修改变量为全局变量
global beta1 beta2 k1 k2 gamma
dx=zeros(1,1);
dx(1)=-gamma*x+beta1*k1*x*(1-x)+beta2*k2*x^2*(1-x);
end
