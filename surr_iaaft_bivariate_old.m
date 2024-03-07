%% Surrogates IAAFT for Bivariate Systems

% DATA: Bivariate Process (N*2)
% iX1= Index First Process 
% iX2= Index Second Process
% nit: Number of Iteractions (default 7)
% stop: if I put 'spe' it comes out with the spectrum preserved, 'dis' it comes out with the distribution preserved

function [xs,ys]=surr_iaaft_bivariate_old(DATA,iX1,iX2,nit,stop)

narginchk(1,5);%min e max of input arguments
if nargin < 5, stop='spe'; end % default (Spectrum its preserved)
if nargin < 4, nit=7; end %default (7 iterations)

x=DATA(:,iX1);
y=DATA(:,iX2);

[xsorted,~]=sort(x);
[ysorted,~]=sort(y);

%% inizializzazione
fx=fft(x);
fy=fft(y);
mx=abs(fx);
my=abs(fy);
fasex_org=angle(fft(x));
fasey_org=angle(fft(y));

p=randperm(size(x,1));
xpermuted=x(p,:);
% ypermuted=y(p,:);

fasex=angle(fft(xpermuted));
% fasey=angle(fft(ypermuted));
phase_diff = fasex - fasex_org;
fasey = fasey_org + phase_diff;

fxs=mx.*(cos(fasex)+1i*sin(fasex));
fys=my.*(cos(fasey)+1i*sin(fasey));
xs=ifft(fxs);xs=real(xs);xs=xs-mean(xs,1);
ys=ifft(fys);ys=real(ys);ys=ys-mean(ys,1);

%[xs, ys, ~, ~]=Shuff_Bivariate(x,y);

%% ciclo

for i=1:nit
    
    % step 1: impone lo spettro
    fasexs=angle(fft(xs));
    fxs=mx.*(cos(fasexs)+1i*sin(fasexs));
    xs=ifft(fxs);
    xs=real(xs);
    xs=xs-mean(xs,1);
    
%     faseys=fasexs;
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys,1);

    % step 2: impone la distribuzione
    [~,xsindice]=sort(xs);
    xpermuted=zeros(length(x),size(x,2));
    for i=1:length(x)
        for ll=1:size(xpermuted,2)
            xpermuted(xsindice(i,ll),ll)=xsorted(i,ll);
        end
    end
    xs=xpermuted;

    
    [~,ysindice]=sort(ys);
    ypermuted=zeros(length(y),size(y,2));
    for i=1:length(y)
        for ll=1:size(ypermuted,2)
            ypermuted(ysindice(i,ll),ll)=ysorted(i,ll);
        end
    end
    ys=ypermuted;

end

%se volevo conservare lo spettro, faccio 1 altro mezzo giro dove impongo solo quello
if strcmp(stop,'spe')
    fasexs=angle(fft(xs));
    fxs=mx.*(cos(fasexs)+1i*sin(fasexs));
    xs=ifft(fxs);xs=real(xs);
    xs=xs-mean(xs,1);
    
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys,1);
end