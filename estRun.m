function [x,y,theta,internalStateOut] = estRun(time, dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
% In this function you implement your estimator. The function arguments
% are:
%  time: current time in [s]
%  dt: current time step [s]
%  internalStateIn: the estimator internal state, definition up to you.
%  steeringAngle: the steering angle of the bike, gamma, [rad]
%  pedalSpeed: the rotational speed of the pedal, omega, [rad/s]
%  measurement: the position measurement valid at the current time step
%
% Note: the measurement is a 2D vector, of x-y position measurement.
%  The measurement sensor may fail to return data, in which case the
%  measurement is given as NaN (not a number).
%
% The function has four outputs:
%  est_x: your current best estimate for the bicycle's x-position
%  est_y: your current best estimate for the bicycle's y-position
%  est_theta: your current best estimate for the bicycle's rotation theta
%  internalState: the estimator's internal state, in a format that can be understood by the next call to this function

% Example code only, you'll want to heavily modify this.
% this needs to correspond to your init function:

%load state variables from internal state object
xhat_prev = internalStateIn.x;
yhat_prev = internalStateIn.y;
theta_hat_prev = internalStateIn.theta;


gamma = steeringAngle;  %steering angle

Pm = internalStateIn.Pm; %variance of states; x,y,xl,yl,theta,gamma

n=3; %number of states

x_var = .05*1; %variance of process noise for x state 
y_var = .05*1; %variance of process noise for y state
theta_var = .1*1; %variance of process noise for theta state

G = diag([x_var,y_var,theta_var]); %variance of process noise 

%measurement noises

p_cov = [1.0893,1.533;1.5333,2.988];


%% Generate 12 sigma points 

B=0.8; % bike baseline

r=.425; % tire radius

v = r*5*pedalSpeed; %calculate bike speed from pedal speed

matrix_root = chol(n*Pm); %find matrix root of Pm so we can make sigma points

state_vec_m = vectorize_state(internalStateIn); %transform internalState into a state vector


%generate sigma points 
for i = 1:n
    s_xm_prev(:,i) = state_vec_m + matrix_root(:,i);
    s_xm_prev(:,i+n) = state_vec_m - matrix_root(:,i);
end



%% Compute prior sigma points 

for i = 1:2*n
    sxp_state(i) = compute_s_points(s_xm_prev(1,i),s_xm_prev(2,i),s_xm_prev(3,i),gamma,dt,B,v);
end
%compute statistics

state_p.x=0; %state_p is the state vector for xhat_p
state_p.y=0;
state_p.theta=0;

%compute averages of sigma points
for i = 1:2*n 
    state_p.x = state_p.x+ sxp_state(i).x/(2*n);
    state_p.y = state_p.y+ sxp_state(i).y/(2*n);
    state_p.theta = state_p.theta+ sxp_state(i).theta/(2*n);
end

%below is the xhat_p vector
state_vec_p = vectorize_state(state_p); %vectorize state object

%compute Pp stats
Pp = zeros(n);

for i = 1:2*n 
    si_vector = vectorize_state(sxp_state(i)); %convert each state to a vector
    Pp = Pp + 1/(2*n)*(si_vector-state_vec_p)*(si_vector-state_vec_p)'; % sum to find Pp   
end

Pp = Pp + G;

%% Calculate measurement update

measurement_available = true; % this tracks whether there is position data avaialbe, the will change the h function

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    px = measurement(1);
    py = measurement(2);
else
    measurement_available = false;
end

%create sz points, expected measurement sz points  


%the measurement update only occurs when we have a valid measurement,
%otherwise it just skips the measurement update and uses the prior stats to
%update the estimator

if measurement_available
    for i=1:2*n
        %compute sz points by pushing sigma points through measurment model
        sz(:,i) = compute_sz_points(sxp_state(i),dt,r,measurement_available,B); %each column is an spoint 
    end

    zhat = sum(sz,2)/(2*n);  %compute average of measurement estimate
    %pzz_forloop = zeros(2);
    
    M=p_cov; %measurement covariance 
    


    for i = 1:2*n
        % removes the mean from sz so we can compute Pzz
        sz_centered(:,i) = sz(:,i) - zhat;
        %pzz_forloop = pzz_forloop + sz_centered(:,i)*sz_centered(:,i)'/(2*n);
    end

    %pzz_forloop = pzz_forloop + M;
    
    Pzz = sz_centered*sz_centered'/(2*n) +M; % compute variance of sz points

    %this for loop removes the mean from prior sigma points
    for i = 1:2*n
        sxp_centered(:,i) = vectorize_state(sxp_state(i))-state_vec_p; %center sxp matrix to remove the mean
    end
    

    Pxz = sxp_centered*sz_centered'/(2*n);    %compute Pxz

    K = Pxz*inv(Pzz); %compute kalman gain

    %compute xm, which is dependent on if we have xy data or not

    z = [measurement(1);measurement(2)];
    xm_vec = state_vec_p + K*(z-zhat);


else
    K=0;
    Pzz=0;
    xm_vec = state_vec_p;
end

%compute Pm

Pm = Pp-K*Pzz*K';


%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:
internalStateIn;

velocity = r*5*pedalSpeed;
finalState = devectorize_state(xm_vec) ;%xhat_m state



internalStateOut.x = finalState.x;
internalStateOut.y = finalState.y;
internalStateOut.theta = finalState.theta;


internalStateOut.Pm = Pm;

x=finalState.x;
y = finalState.y;
theta = finalState.theta;

end

%% Functions 
function [state] = compute_s_points(x,y,theta,gamma,dt,B,v)

%computes model update with s points
state.x = v*cos(theta)*dt + x;    
state.y = v*sin(theta)*dt + y;
state.theta = theta + v/B*tan(gamma)*dt;




end

function [state_vector] = vectorize_state(state)
    %utility function that turns a state object into a state vector
    state_vector(1) = state.x;
    state_vector(2) = state.y;

    state_vector(3) = state.theta;

    state_vector = state_vector';
end

function [meas_vector] = compute_sz_points(state_p_i,Ts,r,meas_available,B,gamma)
    %computes measurement model update with s points
   if meas_available
       
       Px = state_p_i.x + 1/2 * B *cos(state_p_i.theta);
       Py = state_p_i.y + 1/2 * B * sin(state_p_i.theta);
       meas_vector = [Px;Py];
       
   
   end
   
end

function [state_object] = devectorize_state(state_vec)
    %utility that turns a state vector into a state object
    state_object.x = state_vec(1);
    state_object.y = state_vec(2);

    state_object.theta = state_vec(3);


end