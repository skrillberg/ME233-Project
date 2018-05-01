function internalState = estInitialize
% Fill in whatever initialization you'd like here. This function
% generates the internal state of the estimator at time 0. You may do
% whatever you like here, but you must return something that is in the
% format as may be used by your run() function.
%

% we make the interal state a structure, with the first three elements the
% positions x, y; the angle theta; and our favourite colour.

% note that there is *absolutely no prescribed format* for this internal state.
% You can put in it whatever you like. Probably, you'll want to keep the position
% and angle, and probably you'll remove the color.
internalState.x = 0;
internalState.y = 0;
internalState.xl = 0;
internalState.yl = 0;
internalState.theta = 45/180*3.14159;


internalState.dt_prev=0.1;
internalState.gamma = 0;

internalState.Pm = eye(6);


end


