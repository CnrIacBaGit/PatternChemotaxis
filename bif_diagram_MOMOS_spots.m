% Bifurcation diagram in the (beta,q) parameter space
% MOMOS - stripes

clear all
close all
clc

% Parameters
Du = 0.6;
Dv = 0.6;
k1 = 0.4;
k2 = 0.6;
c = 0.8;

% Conditions fore chemotaxis driven instability
c1 = @(beta,q) -k2*Du+Dv*(-k1-2*sqrt(c*q))+beta*k1*sqrt(c/q);
f_pattern = @(beta,q) [c1(beta,q)]^2-8*Dv*Du*k2*sqrt(c*q);

vec_beta = [0:1e-4:2];
vec_q = [0:1e-4:0.1];
m1 = length(vec_beta);
m2 = length(vec_q);
val = zeros(m1,m2);
for i = 1 : length(vec_beta)
    beta = vec_beta(i);
    for j = 1 : length(vec_q)
        q = vec_q(j);
        val1 = c1(beta,q);
        if val1 > 0
            value = f_pattern(beta,q);
            if value > 0
                val(i,j) = 1;
            end
        end
    end
end


[bb,qq] = meshgrid(vec_beta,vec_q);
figure
pcolor(bb,qq,val')
shading interp
xlabel('\beta')
ylabel('q')
hold on
plot(1.2,0.075,'*r')