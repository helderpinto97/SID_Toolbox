%% Surrogates IAAFT for Bivariate Systems

% DATA: Bivariate Process (N*2)
% iX1= Index First Process 
% iX2= Index Second Process
% nit: Number of Iteractions (default 7)
% stop: 'spe' preserves the spectrum, 'dis' preserves the distribution

function [xs,ys]=surr_IaaftBivariate(DATA,iX1,iX2,nit,stop)

narginchk(1,5);
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
fasex_org=angle(fx);
fasey_org=angle(fy);

xpermuted=surr_ShuffColumn(x);  % shuffling

fasex=angle(fft(xpermuted));
phase_diff = fasex - fasex_org;
fasey = fasey_org + phase_diff;

fxs=mx.*(cos(fasex)+1i*sin(fasex));
fys=my.*(cos(fasey)+1i*sin(fasey));
xs=ifft(fxs);xs=real(xs);xs=xs-mean(xs,1);
ys=ifft(fys);ys=real(ys);ys=ys-mean(ys,1);

%% ciclo

for i=1:nit
    
    % step 1: imposing the spectrum
    fasexs=angle(fft(xs));
    fxs=mx.*(cos(fasexs)+1i*sin(fasexs));
    xs=ifft(fxs);
    xs=real(xs);
    xs=xs-mean(xs,1);

    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys,1);

    % step 2: imposing the distribution
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

%imposing the same spectrum
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