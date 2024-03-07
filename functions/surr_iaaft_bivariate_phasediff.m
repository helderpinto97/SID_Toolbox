%% Surrogates IAAFT for Bivariate Systems

% DATA: Bivariate Process (N*2)
% iX1= Index First Process 
% iX2= Index Second Process
% nit: Number of Iteractions (default 7)
% stop: if I put 'spe' it comes out with the spectrum preserved, 'dis' it comes out with the distribution preserved

function [xs,ys]=surr_iaaft_bivariate_phasediff(DATA,iX1,iX2,nit,stop)

%error(nargchk(1,4,nargin));%min e max di input arguments
if nargin < 5, stop='spe'; end %default matcha lo spettro
if nargin < 4, nit=7; end %default 7 iterazioni

% clear;close all;
% percorso='D:\johnny\lavoro\integrate_nlpred\elaborati_loo_si\';% percorso dei dati da analizzare
% nomefile='b-ca.prn';
% rs=load([percorso nomefile]);
% y=rs(:,1);
% y=(y-mean(y))/std(y);
x=DATA(:,[iX1]);
y=DATA(:,[iX2]);

[xsorted,~]=sort(x);
[ysorted,~]=sort(y);

%% inizializzazione
fx=fft(x);
fy=fft(y);
mx=abs(fx);
my=abs(fy);
fasex_org=angle(fx);
fasey_org=angle(fy);

p=randperm(size(x,1));
xs=x(p);

%phase_diff=angle(fft(y(p))) - angle(fft(x(p)));

phase_diff = fasey_org - fasex_org;

%% ciclo
for i=1:nit

    % step 1: impone lo spettro
    fasexs=angle(fft(xs));
    fxs=mx.*exp(1i*fasexs);%(cos(fasexs)+1i*sin(fasexs));
    xs=ifft(fxs);
    xs=real(xs);
    %xs=xs-mean(xs);

    faseys = fasexs + phase_diff;

%     faseys=angle(fft(ys));
    fys=my.*exp(1i*faseys);%(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    
%     faseys=fasexs;
    % ;
    

    
    %ys=ys-mean(ys);

    % step 2: impone la distribuzione
    [~,xsindice]=sort(xs);

%     for ii=1:length(x)
%         xs(xsindice(ii))=xsorted(ii);
%     end

    xs(xsindice) = xsorted;
    
    
    [~,ysindice]=sort(ys);
%     for ii=1:length(y)
%         ys(ysindice(ii))=ysorted(ii);
%     end

    ys(ysindice) = ysorted;

end

%se volevo conservare lo spettro, faccio 1 altro mezzo giro dove impongo solo quello
if stop=='spe'
    fasexs=angle(fft(xs));
    fxs=mx.*exp(1i*fasexs);%(cos(fasexs)+1i*sin(fasexs));
    xs=ifft(fxs);xs=real(xs);
%     xs=xs-mean(xs);
    
    %faseys=angle(fft(ys));
    faseys = fasexs + phase_diff;
    fys=my.*exp(1i*faseys);%(cos(faseys)+1i*sin(faseys));
    ys=ifft(fys);ys=real(ys);
%     ys=ys-mean(ys);
end
