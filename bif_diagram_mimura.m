% Bifurcation diagram in the (beta,q) parameter space
% The Mimura-Tsujikawa model - hexagons

Du = 0.0625;
Dv = 1;
k1 = 1; 
k2 = 32;

% Conditions for chemotaxis driven instability
c1 = @(beta,q) beta*k1-Du*k2-Dv*q;
f_pattern = @(beta,q) [c1(beta,q)]^2-4*Dv*Du*q*k2;

vec_beta = [0:1e-2:20];
vec_q = [0:1e-2:10];
m1 = length(vec_beta);
m2 = length(vec_q);
val = -1*ones(m1,m2);
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
plot(17,7,'*r')
