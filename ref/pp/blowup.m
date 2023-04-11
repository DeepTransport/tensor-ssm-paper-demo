function [value,isterminal,direction] = blowup(~, y)
% the event function for ode solver
% terminate when the states blow up

% value = ~((y(1)>1e3)||(y(2)>1e3)||(y(1)<0)||(y(2)<0));
value = ~((y(1)>5e2)||(y(2)>2e2));
isterminal = 1;
direction = 0;

end