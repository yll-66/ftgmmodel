# ftgmmodel
clc
clear all
load('Data_dep_store_y.mat') %70 dep_store
load('origin_valid') % split origin and valid set
%%
for arr=1:70 %%Cate_paroduct
    data=Data_dep_store_y(:,arr);
    W = ; % window size = origin size
    [k,o]=size(data);
    wid=2:1:k;
    f =; % validation data
    ff=15; % test data length
    m=k-W-f-ff+1;
    t=[1:1:k-W+1]';
    nub=length(t);
    r=1;%时间间隔
    w=2*pi/12;%角转率
    
    %% 生成original data
    XX=[];
    for ii=1:1:k-W+1
         XX(ii,:) = data(ii:ii+W-1);
    end
    X=[];
    for i=1:1:m
         X(i,:) = data(i:i+W-1);
    end
    Y1=[X(1,:);
       cumsum(r.*X(2:m,:))+X(1,:)];
    U=ones(m-1,1);
    %%    %%灰色参数估计%%%%%%%
     for j=1:1:W
          % N=1傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
            pa_1{j}=0.5*(Y1(1:m-1,j)+Y1(2:m,j));
            pa_3{j}=0.5*[Y1(1:m-1,j).*cos(w*t(1:m-1))+Y1(2:m,j).*cos(w*t(2:m)) Y1(1:m-1,j).*sin(w*t(1:m-1))+Y1(2:m,j).*sin(w*t(2:m))];
            pa_4{j}=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m))];
            A{j}=[pa_1{j} U pa_3{j} pa_4{j}];
            pa(:,j)=inv(A{j}'*A{j})*A{j}'*X(2:m,j);
            %% %%%%Origin-based 参数估计%%%%%%%%%%%
            %Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
            Z(:,j)=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*((pa(1,j)+sum(pa(3:4,j)))*X(1)+pa(2,j)+sum(pa(5:6,j)));
            F1=[cos(w*t(2:m)) sin(w*t(2:m))];%%%%%%N=1系数
            S{j}=Z(:,j).*F1;%%%%%%N=1系数
            B{j}=[Z(:,j) U S{j} F1];%矩阵
            P{j}=inv(B{j}'*B{j})*B{j}'*X(2:m,j);%求解的模型参数
            %eta{j}=X(1,j);
            eta(:,j) =(pa(1,j)+sum(pa(3:4,j)))*X(1,j)+pa(2,j)+sum(pa(5:6,j));
            x0{j} = [X(1,j) eta(:,j)];%初始值选择
            %%  %%%%%积分匹配模型拟合误差
            [t,x{j}] = ode45(@(t,x) odefcn(t,x,P{j},w), t, x0{j});
            x_fit(:,j)=x{j}(:,2);
         
           %% N=2傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
%             pa_1{j}=0.5*(Y1(1:m-1,j)+Y1(2:m,j));
%             pa_3{j}=0.5*[Y1(1:m-1,j).*cos(w*t(1:m-1))+Y1(2:m,j).*cos(w*t(2:m)) Y1(1:m-1,j).*sin(w*t(1:m-1))+Y1(2:m,j).*sin(w*t(2:m))];
%             pa_4{j}=0.5*[Y1(1:m-1,j).*cos(2*w*t(1:m-1))+Y1(2:m,j).*cos(2*w*t(2:m)) Y1(1:m-1,j).*sin(2*w*t(1:m-1))+Y1(2:m,j).*sin(2*w*t(2:m))];
%             pa_5{j}=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m)) cos(2*w*t(1:m-1))+cos(2*w*t(2:m)) sin(2*w*t(1:m-1))+sin(2*w*t(2:m))];
%             A{j}=[pa_1{j} U pa_3{j} pa_4{j} pa_5{j}];
%             pa(:,j)=inv(A{j}'*A{j})*A{j}'*X(2:m,j);
%             %% %%%%Origin-based 参数估计%%%%%%%%%%%
%             %Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
%             Z(:,j)=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*((pa(1,j)+sum(pa(3:6,j)))*X(1)+pa(2,j)+sum(pa(7:10,j)));
%             F1=[cos(w*t(2:m)) sin(w*t(2:m)) cos(2*w*t(2:m)) sin(2*w*t(2:m))];%%%%%%N=2系数
%             S{j}=Z(:,j).*F1;%%%%%%N=1系数
%             B{j}=[Z(:,j) U S{j} F1];%矩阵
%             P{j}=inv(B{j}'*B{j})*B{j}'*X(2:m,j);%求解的模型参数
%             %eta{j}=X(1,j);
%             eta(:,j) =(pa(1,j)+sum(pa(3:6,j)))*X(1,j)+pa(2,j)+sum(pa(7:10,j));
%             x0{j} = [X(1,j) eta(:,j)];%初始值选择
%             %%  %%%%%积分匹配模型拟合误差
%             [t,x{j}] = ode45(@(t,x) odefcnn(t,x,P{j},w), t, x0{j});
%             x_fit(:,j)=x{j}(:,2);
%     
           %% N=3傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
%             pa_1{j}=0.5*(Y1(1:m-1,j)+Y1(2:m,j));
%             pa_3{j}=0.5*[Y1(1:m-1,j).*cos(w*t(1:m-1))+Y1(2:m,j).*cos(w*t(2:m)) Y1(1:m-1,j).*sin(w*t(1:m-1))+Y1(2:m,j).*sin(w*t(2:m))];
%             pa_4{j}=0.5*[Y1(1:m-1,j).*cos(2*w*t(1:m-1))+Y1(2:m,j).*cos(2*w*t(2:m)) Y1(1:m-1,j).*sin(2*w*t(1:m-1))+Y1(2:m,j).*sin(2*w*t(2:m))];
%             pa_5{j}=0.5*[Y1(1:m-1,j).*cos(3*w*t(1:m-1))+Y1(2:m,j).*cos(3*w*t(2:m)) Y1(1:m-1,j).*sin(3*w*t(1:m-1))+Y1(2:m,j).*sin(3*w*t(2:m))];
%             pa_6{j}=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m)) cos(2*w*t(1:m-1))+cos(2*w*t(2:m)) sin(2*w*t(1:m-1))+sin(2*w*t(2:m)) cos(3*w*t(1:m-1))+cos(3*w*t(2:m)) sin(3*w*t(1:m-1))+sin(3*w*t(2:m))];
%             A{j}=[pa_1{j} U pa_3{j} pa_4{j} pa_5{j} pa_6{j}];
%             pa(:,j)=inv(A{j}'*A{j})*A{j}'*X(2:m,j);
%             %% %%%%Origin-based 参数估计%%%%%%%%%%%
%             %Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
%             Z(:,j)=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*((pa(1,j)+sum(pa(3:8,j)))*X(1)+pa(2,j)+sum(pa(9:14,j)));
%             F1=[cos(w*t(2:m)) sin(w*t(2:m)) cos(2*w*t(2:m)) sin(2*w*t(2:m)) cos(3*w*t(2:m)) sin(3*w*t(2:m))];%%%%%%N=3系数
%             S{j}=Z(:,j).*F1;%%%%%%N=1系数
%             B{j}=[Z(:,j) U S{j} F1];%矩阵
%             P{j}=inv(B{j}'*B{j})*B{j}'*X(2:m,j);%求解的模型参数
%             %eta{j}=X(1,j);
%             eta(:,j) =(pa(1,j)+sum(pa(3:8,j)))*X(1,j)+pa(2,j)+sum(pa(9:14,j));
%             x0{j} = [X(1,j) eta(:,j)];%初始值选择
%             %%  %%%%%积分匹配模型拟合误差
%             [t,x{j}] = ode45(@(t,x) odefcnnn(t,x,P{j},w), t, x0{j});
%             x_fit(:,j)=x{j}(:,2);
%     
            %% N=4傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
%             pa_1{j}=0.5*(Y1(1:m-1,j)+Y1(2:m,j));
%             pa_3{j}=0.5*[Y1(1:m-1,j).*cos(w*t(1:m-1))+Y1(2:m,j).*cos(w*t(2:m)) Y1(1:m-1,j).*sin(w*t(1:m-1))+Y1(2:m,j).*sin(w*t(2:m))];
%             pa_4{j}=0.5*[Y1(1:m-1,j).*cos(2*w*t(1:m-1))+Y1(2:m,j).*cos(2*w*t(2:m)) Y1(1:m-1,j).*sin(2*w*t(1:m-1))+Y1(2:m,j).*sin(2*w*t(2:m))];
%             pa_5{j}=0.5*[Y1(1:m-1,j).*cos(3*w*t(1:m-1))+Y1(2:m,j).*cos(3*w*t(2:m)) Y1(1:m-1,j).*sin(3*w*t(1:m-1))+Y1(2:m,j).*sin(3*w*t(2:m))];
%             pa_6{j}=0.5*[Y1(1:m-1,j).*cos(4*w*t(1:m-1))+Y1(2:m,j).*cos(4*w*t(2:m)) Y1(1:m-1,j).*sin(4*w*t(1:m-1))+Y1(2:m,j).*sin(4*w*t(2:m))];
%             pa_7{j}=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m)) cos(2*w*t(1:m-1))+cos(2*w*t(2:m)) sin(2*w*t(1:m-1))+sin(2*w*t(2:m)) cos(3*w*t(1:m-1))+cos(3*w*t(2:m)) sin(3*w*t(1:m-1))+sin(3*w*t(2:m)) cos(4*w*t(1:m-1))+cos(4*w*t(2:m)) sin(4*w*t(1:m-1))+sin(4*w*t(2:m))];
%             A{j}=[pa_1{j} U pa_3{j} pa_4{j} pa_5{j} pa_6{j} pa_7{j}];
%             pa(:,j)=inv(A{j}'*A{j})*A{j}'*X(2:m,j);
%             %% %%%%Origin-based 参数估计%%%%%%%%%%%
%             %Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
%             Z(:,j)=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*((pa(1,j)+sum(pa(3:10,j)))*X(1)+pa(2,j)+sum(pa(11:18,j)));
%             F1=[cos(w*t(2:m)) sin(w*t(2:m)) cos(2*w*t(2:m)) sin(2*w*t(2:m)) cos(3*w*t(2:m)) sin(3*w*t(2:m)) cos(4*w*t(2:m)) sin(4*w*t(2:m))];%%%%%%N=4系数
%             S{j}=Z(:,j).*F1;%%%%%%N=1系数
%             B{j}=[Z(:,j) U S{j} F1];%矩阵
%             P{j}=inv(B{j}'*B{j})*B{j}'*X(2:m,j);%求解的模型参数
%             %eta{j}=X(1,j);
%             eta(:,j) =(pa(1,j)+sum(pa(3:10,j)))*X(1,j)+pa(2,j)+sum(pa(11:18,j));
%             x0{j} = [X(1,j) eta(:,j)];%初始值选择
%             %%  %%%%%积分匹配模型拟合误差
%             [t,x{j}] = ode45(@(t,x) odefcnnnn(t,x,P{j},w), t, x0{j});
%             x_fit(:,j)=x{j}(:,2);
%     
      %% N=5傅里叶级数%%%%%%%%%Origin-based 参数估计%%%%%%%%%%%
%             pa_1{j}=0.5*(Y1(1:m-1,j)+Y1(2:m,j));
%             pa_3{j}=0.5*[Y1(1:m-1,j).*cos(w*t(1:m-1))+Y1(2:m,j).*cos(w*t(2:m)) Y1(1:m-1,j).*sin(w*t(1:m-1))+Y1(2:m,j).*sin(w*t(2:m))];
%             pa_4{j}=0.5*[Y1(1:m-1,j).*cos(2*w*t(1:m-1))+Y1(2:m,j).*cos(2*w*t(2:m)) Y1(1:m-1,j).*sin(2*w*t(1:m-1))+Y1(2:m,j).*sin(2*w*t(2:m))];
%             pa_5{j}=0.5*[Y1(1:m-1,j).*cos(3*w*t(1:m-1))+Y1(2:m,j).*cos(3*w*t(2:m)) Y1(1:m-1,j).*sin(3*w*t(1:m-1))+Y1(2:m,j).*sin(3*w*t(2:m))];
%             pa_6{j}=0.5*[Y1(1:m-1,j).*cos(4*w*t(1:m-1))+Y1(2:m,j).*cos(4*w*t(2:m)) Y1(1:m-1,j).*sin(4*w*t(1:m-1))+Y1(2:m,j).*sin(4*w*t(2:m))];
%             pa_7{j}=0.5*[Y1(1:m-1,j).*cos(5*w*t(1:m-1))+Y1(2:m,j).*cos(5*w*t(2:m)) Y1(1:m-1,j).*sin(5*w*t(1:m-1))+Y1(2:m,j).*sin(5*w*t(2:m))];
%             pa_8{j}=0.5*[cos(w*t(1:m-1))+cos(w*t(2:m)) sin(w*t(1:m-1))+sin(w*t(2:m)) cos(2*w*t(1:m-1))+cos(2*w*t(2:m)) sin(2*w*t(1:m-1))+sin(2*w*t(2:m)) cos(3*w*t(1:m-1))+cos(3*w*t(2:m)) sin(3*w*t(1:m-1))+sin(3*w*t(2:m)) cos(4*w*t(1:m-1))+cos(4*w*t(2:m)) sin(4*w*t(1:m-1))+sin(4*w*t(2:m)) cos(5*w*t(1:m-1))+cos(5*w*t(2:m)) sin(5*w*t(1:m-1))+sin(5*w*t(2:m))];
%             A{j}=[pa_1{j} U pa_3{j} pa_4{j} pa_5{j} pa_6{j} pa_7{j} pa_8{j}];
%             pa(:,j)=inv(A{j}'*A{j})*A{j}'*X(2:m,j);
%             %% %%%%Origin-based 参数估计%%%%%%%%%%%
%             %Z{j}=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*X(1,j);
%             Z(:,j)=0.5*Y1(1:m-1,j)+0.5*Y1(2:m,j)+0.5*r*((pa(1,j)+sum(pa(3:12,j)))*X(1)+pa(2,j)+sum(pa(13:22,j)));
%             F1=[cos(w*t(2:m)) sin(w*t(2:m)) cos(2*w*t(2:m)) sin(2*w*t(2:m)) cos(3*w*t(2:m)) sin(3*w*t(2:m)) cos(4*w*t(2:m)) sin(4*w*t(2:m)) cos(5*w*t(2:m)) sin(5*w*t(2:m))];%%%%%%N=5系数
%             S{j}=Z(:,j).*F1;%%%%%%N=5系数
%             B{j}=[Z(:,j) U S{j} F1];%矩阵
%             P{j}=inv(B{j}'*B{j})*B{j}'*X(2:m,j);%求解的模型参数
%             %eta{j}=X(1,j);
%             eta(:,j) =(pa(1,j)+sum(pa(3:12,j)))*X(1,j)+pa(2,j)+sum(pa(13:22,j));
%             x0{j} = [X(1,j) eta(:,j)];%初始值选择
%             %%  %%%%%积分匹配模型拟合误差
%             [t,x{j}] = ode45(@(t,x) odefcnnnnn(t,x,P{j},w), t, x0{j});
%             x_fit(:,j)=x{j}(:,2);
%     
     end
    
     %% accuracy metrics
    for l=1:W
        %%%%fitting
        error_fit_1(:,l)=X(2:m,l)-XX(1:m-1,l);
%         error_fit_2(:,l)=XX(2:m+f,l)-XX(1:m+f-1,l);
        error_fit_2(:,l)=XX(2:m+f,l)-x_fit(1:m+f-1,l);
        %%%%fitted 
        error_fit(:,l)=X(1:m,l)-x_fit(1:m,l);
        pe_fit(:,l)= (error_fit(:,l))./X(1:m,l)*100;
        %%%%%alidation 
        error_valid(:,l)=XX(m+1:m+f,l)-x_fit(m+1:m+f,l);
        pe_valid(:,l)= (error_valid(:,l))./XX(m+1:m+f,l)*100;
        %%%%%prediction
        error_fore(:,l)=XX(m+f+1:m+f+ff,l)-x_fit(m+f+1:m+f+ff,l);
        pe_fore(:,l)= (error_fore(:,l))./XX(m+f+1:m+f+ff,l)*100;
        spe_fore(:,l)= 2*(error_fore(:,l))./(abs(XX(m+f+1:m+f+ff,l))+abs(x_fit(m+f+1:m+f+ff,l)))*100;
        RMSE_fit_w(:,l) = sqrt(mean(mean((error_fit(:,l)).^2,2)));
        RMSE_valid_w(:,l) = sqrt(mean(mean((error_valid(:,l)).^2,2)));
    end

    a_RMSE_vad(arr,:) = [RMSE_valid_w]';
    a_RMSE_fit(arr,:) = [RMSE_fit_w]';
    
    %% %%%each-step RMSE%%%%%validation
    RMSE_fit_step = sqrt(mean((error_fit).^2,2));
    RMSE_fit_sum = sqrt(mean(mean((error_fit).^2,2)));
    RMSE_valid_step = sqrt(mean((error_valid).^2,2));
    RMSE_valid_sum = sqrt(mean(mean((error_valid).^2,2)));
    RMSE_fore_step = sqrt(mean((error_fore).^2,2));
    RMSE_fore_sum = sqrt(mean(mean((error_fore).^2,2)));

    RMSSE_fore_step = sqrt(mean((error_fore).^2 ./mean((error_fit_2).^2,1),2));
    RMSSE_fore_sum = sqrt(mean(mean((error_fore).^2 ./mean((error_fit_2).^2,1),2)));

    MASE_fore_step = mean(abs(error_fore)./mean(abs(error_fit_2),1),2);
    MASE_fore_sum = mean(MASE_fore_step);

    error_ff(arr,:) = [RMSE_fit_sum MASE_fit_sum RMSSE_fit_sum]';
    
    error_f(arr,:) = [RMSE_fore_sum MASE_fore_sum RMSSE_fore_sum]';
    
    RMSE_f_step(arr,:) = [RMSE_fore_step'];
    MASE_f_step(arr,:) = [MASE_fore_step'];
    RMSSE_f_step(arr,:) = [RMSSE_fore_step'];
    e_f_step(arr,:) = [RMSE_f_step(arr,:) zeros(1,1)  MASE_f_step(arr,:) zeros(1,1) RMSSE_f_step(arr,:)];

end
