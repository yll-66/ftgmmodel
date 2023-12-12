function [x_f,rmseX,rmseX_1_step] = ftgm_1(ts,xx,xs,xf,nf,m,r,w)

% ts: time instant series 
% xx:  original time series
% xs: training data
% m: fitting steps
% nf: forcasting steps
% r: time interval in ts 
% w: angular rotation rate

% P: estimated parameters, not include 0 item
% x_fit: fitted series corresponding to x
% x_fore:   forecasting series with length nf 

%% N=1 Fourier order

    % cumulative sum
    Y1 = cumsum([1;diff(ts(1:end-nf))].*xx); 

    U = ones(m-1,1);
    %%      %%%%%%%%%Origin-based grey model%%%%%%%%%%%

    % close neighbour means for input matrix
    pa_1 = 0.5*(Y1(1:m-1)+Y1(2:m));
    pa_3 = 0.5*[Y1(1:m-1).*cos(w*ts(1:m-1))+Y1(2:m).*cos(w*ts(2:m)) Y1(1:m-1).*sin(w*ts(1:m-1))+Y1(2:m).*sin(w*ts(2:m))];
    pa_4 = 0.5*[cos(w*ts(1:m-1))+cos(w*ts(2:m)) sin(w*ts(1:m-1))+sin(w*ts(2:m))];
    A = [pa_1 U pa_3 pa_4];

    % least squares estimation
    pa = inv(A'*A)*A'*xx(2:m);

    %%       %%%%intergal matching-based model %%%%%%%%%%%

    % Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
    Z = 0.5*Y1(1:m-1)+0.5*Y1(2:m)+0.5*r*((pa(1)+sum(pa(3:4)))*xx(1)+pa(2)+sum(pa(5:6)));
    F1 = [cos(w*ts(2:m)) sin(w*ts(2:m))];
    S = Z.*F1;
    B = [Z U S F1];
    
    % estimated parameters
    P = inv(B'*B)*B'*xx(2:m);

    % initial value for intergal matching-based model 
    %eta{j}=X(1,j);
    eta = (pa(1)+sum(pa(3:4)))*xx(1)+pa(2)+sum(pa(5:6));

    % initial values for Origin-based grey model and intergal matching-based model 
    X0 = [xx(1) eta];

    % fitting and forecasting results
    [t,X] = ode45(@(t,X) odefcn(t,X,P,w), ts, X0);

    % Series fitting
    x_f = X(:,2);
    x_fit = x_f(1:m);

    % Series forecasting
    x_fore = x_f(m+1:end);

    % RMSE error computing
    % fitting errors
    rmseIn = sqrt(sum((xs-x_fit).^2)/m); 

    % forecasting errors
    rmseOut = sqrt(sum((xf-x_fore).^2)/12); 

    rmseX = [rmseIn rmseOut];















