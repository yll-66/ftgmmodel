function [x_f,rmseX,rmseX_5_step] = ftgm_5(t,xx,xs,xf,nf,m,r,w)

% ts: time instant series 
% x:  original time series
% nf: forcasting steps
% h: time interval in ts 

% Pi: estimated parameters, not include 0 item
% xhat: fitted series corresponding to x
% xf:   forecasting series with length nf 

            Y1 = cumsum([1;diff(t(1:end-nf))].*xx); % time iterval
            U=ones(m-1,1);
            %%N=4傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
            %% N=5傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
            pa_1=0.5*(Y1(1:m-1)+Y1(2:m));
            pa_3=0.5*[Y1(1:m-1).*cos(w*t(1:m-1))+Y1(2:m).*cos(w*t(2:m)) Y1(1:m-1).*sin(w*t(1:m-1))+Y1(2:m).*sin(w*t(2:m))];
            pa_4=0.5*[Y1(1:m-1).*cos(2*w*t(1:m-1))+Y1(2:m).*cos(2*w*t(2:m)) Y1(1:m-1).*sin(2*w*t(1:m-1))+Y1(2:m).*sin(2*w*t(2:m))];
            pa_5=0.5*[Y1(1:m-1).*cos(3*w*t(1:m-1))+Y1(2:m).*cos(3*w*t(2:m)) Y1(1:m-1).*sin(3*w*t(1:m-1))+Y1(2:m).*sin(3*w*t(2:m))];
            pa_6=0.5*[Y1(1:m-1).*cos(4*w*t(1:m-1))+Y1(2:m).*cos(4*w*t(2:m)) Y1(1:m-1).*sin(4*w*t(1:m-1))+Y1(2:m).*sin(4*w*t(2:m))];
            pa_7=0.5*[Y1(1:m-1).*cos(5*w*t(1:m-1))+Y1(2:m).*cos(5*w*t(2:m)) Y1(1:m-1).*sin(5*w*t(1:m-1))+Y1(2:m).*sin(5*w*t(2:m))];
            pa_8=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m)) cos(2*w*t(1:m-1))+cos(2*w*t(2:m)) sin(2*w*t(1:m-1))+sin(2*w*t(2:m)) cos(3*w*t(1:m-1))+cos(3*w*t(2:m)) sin(3*w*t(1:m-1))+sin(3*w*t(2:m)) cos(4*w*t(1:m-1))+cos(4*w*t(2:m)) sin(4*w*t(1:m-1))+sin(4*w*t(2:m)) cos(5*w*t(1:m-1))+cos(5*w*t(2:m)) sin(5*w*t(1:m-1))+sin(5*w*t(2:m))];
            A=[pa_1 U pa_3 pa_4 pa_5 pa_6 pa_7 pa_8];
            pa(:)=inv(A'*A)*A'*xx(2:m);
            %% %%%%Origin-based 参数估计%%%%%%%%%%%
            %Z=0.5*Y1(1:m-1)+0.5*Y1(2:m)+0.5*r*X(1);
            Z(:)=0.5*Y1(1:m-1)+0.5*Y1(2:m)+0.5*r*((pa(1)+sum(pa(3:12)))*xx(1)+pa(2)+sum(pa(13:22)));
            F1=[cos(w*t(2:m)) sin(w*t(2:m)) cos(2*w*t(2:m)) sin(2*w*t(2:m)) cos(3*w*t(2:m)) sin(3*w*t(2:m)) cos(4*w*t(2:m)) sin(4*w*t(2:m)) cos(5*w*t(2:m)) sin(5*w*t(2:m))];%%%%%%N=5系数
            S=Z(:).*F1;%%%%%%N=5系数
            B=[Z(:) U S F1];%矩阵
            P=inv(B'*B)*B'*xx(2:m);%求解的模型参数
            eta(:) =(pa(1)+sum(pa(3:12)))*xx(1)+pa(2)+sum(pa(13:22));
            X0 = [xx(1) eta(:)];%初始值选择
            %%  %%%%%积分匹配模型拟合误差
            [t,X] = ode45(@(t,X) odefcnnnnn(t,X,P,w), t, X0);
            x_f=X(:,2);
            x_fit=x_f(1:m);
            x_fore=x_f(m+1:end);
            rmseIn = sqrt(sum((xs-x_fit).^2)/m);               % fitting errors
            rmseOut = sqrt(sum((xf-x_fore).^2)/12); 
            rmseOut_1 = sqrt(sum((xf(end-12+1)-x_fore(end-12+1)).^2)/1); % forecasting errors
            rmseOut_6 = sqrt(sum((xf(end-6+1)-x_fore(end-6)).^2)/1);
            rmseOut_12 = sqrt(sum((xf(end)-x_fore(end)).^2)/1); 
            rmseX = [rmseIn rmseOut];
            rmseX_5_step = [rmseIn rmseOut_1 rmseOut_6 rmseOut_12 rmseOut];













