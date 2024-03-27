function x=Hennon_Map_AdNoise(N,N_transient,delta)

% Delta denote the variance of noise

% Array to store results
x=nan(1,N+N_transient);

% Random Initial Condition
% x(1:2)=randn(1,2);

% Map Parameters
a=1.4;
b=0.3;

flag=true;
nan_counter=0;

while flag

    % Random Initial Condition
    x(1:2)=rand(2,1);

    % Loop
    for i=3:(N_transient+N)
        x(i)=1-a*(x(i-1)^2)+b*x(i-2);
    end

    % Remove transient data
    x=x(N_transient+1:end);
    x=x+delta.*randn(1,N);

    % Check if the generated series contain INF and NAN values
    if ~sum(isnan(x)) && ~sum(isinf(x))
        flag=false;
    elseif nan_counter>5000
        error('Non Stable Dynamics.')
    else
        nan_counter=nan_counter+1;
    end
end
x=x';

end