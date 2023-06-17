function [sigm] = otus(image)
I = image; 
[h,w]=size(I);
xs=zeros(1,256);
kt=10;
for i=1:h
    for j=1:w
       k=I(i,j);
       xs(k+1)=xs(k+1)+1;                
    end
end
%%
pi=zeros(1,255);
mG=0;
%mG 全局的灰度值
 for i = 1:255
        pi(i)=xs(i)/h/w;
        mG=mG+i*pi(i);
 end
 
sigmb=zeros(1,255);
for kt = 1:255
    mk=0;p1=0;
 %P1(k)
 for i = 1:255
     if(i<=kt)
        p1 = p1+pi(i);
     end
 end
 %求m(k)
         for i = 1:kt
             mk = mk+i*pi(i);
         end
p2=1-p1;
sigmb(kt)=((mG*p1-mk)^2)/p1/p2;
end
%%
[m,sigm]=max(sigmb);


end

