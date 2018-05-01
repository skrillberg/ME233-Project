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
xl_hat_prev = internalStateIn.xl;
yl_hat_prev = internalStateIn.yl;
dt_prev = internalStateIn.dt_prev;
gamma_prev = internalStateIn.gamma;

Pm = internalStateIn.Pm; %variance of states; x,y,xl,yl,theta,gamma

n=6;

G = eye(6); %variance of process noise 

%generate 12 sigma points 

B=1.3;
r=.34;
matrix_root = chol(n*Pm);

state_vec_m = vectorize_state(internalStateIn)



for i = 1:n
    s_xm_prev(:,i) = state_vec_m + matrix_root(:,i);
    s_xm_prev(:,i+n) = state_vec_m - matrix_root(:,i);
end

display(s_xm_prev)

%compute prior points 

for i = 1:2*n
    sxp_state(i) = compute_s_points(s_xm_prev(1,i),s_xm_prev(2,i),s_xm_prev(3,i),s_xm_prev(4,i),s_xm_prev(5,i),s_xm_prev(6,i),dt,dt,B);
end
%compute statistics

state_p.x=0; %state_p is the state vector for xhatp
state_p.y=0;
state_p.yl=0;
state_p.xl=0;
state_p.theta=0;
state_p.gamma = 0;
for i = 1:2*n 
    state_p.x = state_p.x+ sxp_state(i).x/(2*n);
    state_p.y = state_p.y+ sxp_state(i).y/(2*n);
    state_p.xl = state_p.xl+ sxp_state(i).xl/(2*n);
    state_p.yl = state_p.yl+ sxp_state(i).yl/(2*n);
    state_p.theta = state_p.theta+ sxp_state(i).theta/(2*n);
    state_p.gamma = state_p.gamma + sxp_state(i).gamma/(2*n);
    
end


state_vec_p = vectorize_state(state_p) %vectorize state vector

%compute Pp stats
Pp = zeros(n);

for i = 1:2*n 
    si_vector = vectorize_state(sxp_state(i));
    Pp = Pp + 1/(2*n)*(si_vector-state_vec_p)*(si_vector-state_vec_p)';    
end

Pp = Pp + G

%calculate measurement update

measurement_available = true; % this tracks whether there is position data avaialbe, the will change the h function

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    px = measurement(1);
    py = measurement(2);
else
    measurement_available = false;
end

%create sz points, expected measurement sz points  

for i=1:2*n
    sz(:,i) = compute_sz_points(sxp_state(i),dt,r,measurement_available,B); %each column is an spoint 
end

zhat = sum(sz,2)/(2*n)  %compute average of measurement estimate

if measurement_available
    pzz_forloop = zeros(4);
else 
    pzz_forloop = zeros(2);
end

for i = 1:2*n
    sz_centered(:,i) = sz(:,i) - zhat;
    pzz_forloop = pzz_forloop + sz_centered(:,i)*sz_centered(:,i)'/(2*n);
end

pzz_forloop
Pzz = sz_centered*sz_centered'/(2*n) % compute variance of sz points

%



%% Functions 
function [state] = compute_s_points(x,y,xl,yl,theta,gamma,dt,dtp,B)
v = sqrt((x  - xl)^2+(y-yl)^2)/dtp;   %calculate velocity 

state.x = v*cos(theta)*dt + x;    
state.y = v*sin(theta)*dt + y;
state.theta = theta + v/B*tan(gamma)*dt;
state.gamma = gamma;
state.xl = x;
state.yl = y;


end

function [state_vector] = vectorize_state(state)
    state_vector(1) = state.x;
    state_vector(2) = state.y;
    state_vector(3) = state.xl;
    state_vector(4) = state.yl;
    state_vector(5) = state.theta;
    state_vector(6) = state.gamma;
    state_vector = state_vector';
end

function [meas_vector] = compute_sz_points(state_p_i,Ts,r,meas_available,B)
    
   if meas_available
       
       gamma_z = state_p_i.gamma;
       omega_z = sqrt( (state_p_i.x-state_p_i.xl)^2 + (state_p_i.y - state_p_i.yl)^2 )/ (5*Ts*r);
       Px = state_p_i.x + 1/2 * B *cos(state_p_i.theta);
       Py = state_p_i.y + 1/2 * B * sin(state_p_i.theta);
       meas_vector = [gamma_z;omega_z;Px;Py]
       
   else
       
       omega_z = sqrt( (state_p_i.x-state_p_i.xl)^2 + (state_p_i.y - state_p_i.yl)^2 )/ (5*Ts*r);
       gamma_z = state_p_i.gamma;
       meas_vector = [gamma_z;omega_z]
       
   end
   
end
%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = x;
internalStateOut.y = y;
internalStateOut.theta = theta;
internalStateOut.color = color;

end