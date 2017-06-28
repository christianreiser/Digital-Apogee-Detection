clc
close all
clear

sigma_baro= 0.001; %35.8229noise baro sensor;
sigma_acc= 0.002;%0.0002 %0.0012 noise acc sensor; 

%% IMU input
fileID = fopen('LOGGER05.TXT','r');
formatSpec = '%s';
Datalogs = fscanf(fileID,formatSpec);
n_datalogs=size(Datalogs);
n_datalogs=n_datalogs(1,2);
n_timesteps= 2390;

%Baro
n=1;
h_baro=zeros(n_timesteps,1);
for i=1:1:n_datalogs %numberof characters in logfile
    if Datalogs(1,i)=='A' && Datalogs(1,i+1)=='l'
        h_baro(n,1)=str2double(Datalogs(1,i+4:i+9));
        n=n+1;
    end
end

%Accelerometer ZK
n=1;
accZK=zeros(n_timesteps,1); 
for i=1:1:n_datalogs %number of characters in logfile
    if Datalogs(1,i)=='c' && Datalogs(1,i+1)=='Z'
        m=3; %measurement has 3 characters: i.e.:9.37
        if Datalogs(1,i+7)~='A'
            m=4; %measurement has 4 characters: i.e.:-9.37
            if Datalogs(1,i+8)~='A'
            m=5; %measurement has 5 characters: i.e.:-10.37
            end
        end
        accZK(n,1)=str2double(Datalogs(1,i+3:1:i+3+m));
        n=n+1;
    end
end

if size(accZK)~=size(h_baro)
	disp('size accZK and h_baro do not match, remove incomplete timestep from logfile')
end

%Magnetometer
alpha= 0; %from magnetometer
beta = 0; %from magnetometer


dt=0.055; %time between measurements in sec
%A is F; m is z; klk-1 is _n;k-1lk-1 is blanc

%% baro mean
first_baro_mean=mean(h_baro(1:10,1)); %mean of first 10 baro measurements
last_baro_mean=mean(h_baro(n_timesteps-10:n_timesteps,1)); %mean of last 10 baro measurements
baro_mean=(first_baro_mean+last_baro_mean)/2; %difference only due to weather?
max(h_baro)
min(h_baro)

%% From accelerated system to inertial system
a_acc=accZK * cos(alpha)* cos(beta); %not finished, add other axes
max_a_acc=max(a_acc);
min_a_acc=min(a_acc);
a_acc_mean=mean(a_acc); 

% mean acceleration of flight sections for plot
takeoff=230;    %right bevore takeoff
asc=625;        %right bevore ascent
apo=1042;       %right bevore apogee
desc=1100;      %right bevore desc
land=1668;      %right bevore landing
last=1897;      %right after landing

acc7=zeros(n_timesteps,1);
acc7(1:takeoff-1,1)=mean(a_acc(1:takeoff-1,1)); %mean of first 10 baro measurements
acc7(takeoff:asc-1,1)=mean(a_acc(takeoff:asc-1,1));
acc7(asc:apo-1,1)=mean(a_acc(asc:apo-1,1));
acc7(apo:desc-1,1)=mean(a_acc(apo:desc-1,1));
acc7(desc:land-1,1)=mean(a_acc(desc:land-1,1));
acc7(land:last-1,1)=mean(a_acc(land:last-1,1));
acc7(last:n_timesteps,1)=mean(a_acc(last:n_timesteps,1));

%% Initialize
z = [h_baro-baro_mean -(a_acc-a_acc_mean)]; %compleatly useless Measurement size(z) %-9.81~a_acc_mean-----------------------change----------------
a=1;
v=1;
h=1;
x=[h;v;a]; %mean state vector
I = eye(3);
P = eye(3);
H= [1 0 0; 0 0 1];
F = [1 dt 1/2 * dt^2; 0 1 dt; 0 0 1]; % maps previous state to next state

%% Measurement noise covariance matrix
% Autocovariance Least-Squares
% http://jbrwww.che.wisc.edu/software/als/

R= [sigma_baro^2 0; 0 sigma_acc^2]; % measurement noise covariance
Q = [0 0 0; 0 0 0; 0 0 1]; % process noise covariance matrix

%% Kalman Gain iteration
for i = 1:20
K = P*H'/(H*P*H' + R); % Kalman gains
P = (eye(3) - K *H)*P;
P = F*P*F' + Q;
end
%display(K)
% display(H)
%display(P)

%% Kalmanfilter
state=zeros(2*n_timesteps,4);
for i=1:n_timesteps
    % Predict (marked as "1" in state)
    x_n = F*x; %size(x_n) %+ B*u; %B*u:updated position due acc; B:control input matrix; u:control vector.
    state(i*2-1,1)=1;
    state(i*2-1,2:1:4)=x_n';
    P_n=F*P*F'+Q; 
    if norm(x_n(2,1))<0.5 && x_n(1,1)>100
        %disp('APOGEE!!')
        %disp(i)
    end
    % Update (marked as "0" in state)
    S= H*P_n*H'+R;
    K=P_n*H'/S;

    y= z(i,:)'-H*x_n;	%size(y)    %trasposed because of array
    x=x_n+K*y;          %size(x_n)
    state(i*2,1)=0;
    state(i*2,2:1:4)=x';
    if norm(x(2,1))<0.5 && x(1,1)>100
        %disp('APOGEE!!')
        %disp(i)
        %disp(h_baro(i,:))
        %disp(accZK(i,:))
    end
    P=(I-K*H)*P_n;
end

figure;
h=h_baro-baro_mean;
t=1:n_timesteps;
plot(t,h,t,state(1:2:2*n_timesteps,2),t,acc7);
ylabel('hight in [m] &acc');
xlabel('timesteps');
legend('barometer','EKF','acc7');

% for i=1:n_timesteps/10
%     m_acc=zeros(n_timesteps,1)
%     a_acc([i*10-10:i+10])=

% %display(state); %1=prediction/0=update;hight;velocity;acceleration
