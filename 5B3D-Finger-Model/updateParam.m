%% 5B3D-Finger-Model: Main
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
% The updateParam function outputs various additinal parameters that are 
% dependent on finger configuration/state:
%   p1 -> contact point 1 along prot_1_imal link aligned with link at_1_is 
%   p2 -> contact point 2 along distal link aligned with link at_1_is 
%   pC -> global position of point connecting links b and c 
%   g -> length of virtual link g
%   tG -> global angle of virtual link g 
%   h -> length of virtual link h
%   tH -> global angle of virtual link g
%   c -> length of compression spring/link c 
%   tC -> global angle of compression spring/link c 
%   n1 -> unit vector for the contact at point 1 (towards link)
%   n2 -> unit vector for the contact at point 2 (towards link)
%   nC-> unit vector for the link C (from O4->O3)

%% This function takes in the following inputs:
%   param -> the structure of parameters to update
%   t_1_ -> angle \theta_1 (scalar)
%   t_2_ -> angle \theta_2 (scalar)
%   t_3_ -> angle \theta_3 (scalar)


function param = updateParam(param, t_1_, t_2_, t_3_)
    param = param;
    % [t_1_, t_2_] = meshgrid(t_1, t_2); % (r, c) = (t_2, t_1)
    
    % These fitted functions come from the ../elipseSolver script
    % p1 = 14.4501062925132 + 28.6147335594148*t_1_ + 9.54112998810529*t_2_ + ...
    %     -14.3693994611985*t_1_.*t_2_ + 2.7779706576232*t_1_.^2 + ...
    %     -2.91159799728414*t_2_.^2;
    % p1 = p1/1000; % convert to m
    % 
    % p2 = 70.5221095238136 + -17.4187426901454*t_1_ + -28.6747238947982*t_2_ + ...
    %     16.0445700081864*t_1_.*t_2_ + -17.171428241453*t_1_.^2.*t_2_ + ...
    %     19.7847306722473*t_1_.*t_2_.^2;
    % p2 = p2/1000; % convert to m
    [p1, p2, a, b] = contactSolver(param, t_1_, t_2_);
    
    pC = param.a*[cos(t_1_); sin(t_1_)] + ...
        param.b*[cos(t_1_+t_2_-param.alpha); sin(t_1_+t_2_-param.alpha)];
    
    % g is onlt_2_ dependent on actuator angle (t_3_ is assumed to be a scalar)
    g = sqrt(param.e^2 + param.d^2 - 2*param.e*param.d*cos(pi-param.gamma-t_3_)); 
    tG = wrapTo2Pi(asin(param.d/g*sin(pi-param.gamma-t_3_)) - param.gamma);
    
    h = sqrt(param.a^2 + param.b^2 - 2*param.a*param.b*cos(pi-param.alpha+t_2_));
    tH = wrapTo2Pi(t_1_ - asin(param.b./h.*sin(pi-param.alpha+t_2_)));
    
    c = sqrt(g.^2 + h.^2 - 2*g.*h.*cos(tH-tG));
    tC = pi + tG - asin(h./c.*sin(tH-tG));

    % contact force unit vectors
    n1 = [cos(t_1_-pi/2); sin(t_1_-pi/2)];
    n2 = [cos(t_1_+t_2_-pi/2); sin(t_1_+t_2_-pi/2)];
    nC = [cos(tC); sin(tC)];

    param.t_1_ = t_1_;
    param.t_2_ = t_2_;
    param.p1 = p1;
    param.p2 = p2;
    param.A = a; % ellipse axis lengths
    param.B = b;
    param.pC = pC;
    param.g = g;
    param.tG = tG;
    param.h = h;
    param.tH = tH;
    param.c = c;
    param.tC = tC;
    param.n1 = n1;
    param.n2 = n2;
    param.nC = nC;
end