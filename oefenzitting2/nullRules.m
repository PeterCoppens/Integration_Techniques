function [] = nullRules(u,f)
% Pas testfuncties 
x = linspace(0,1, length(u));
e(k) = u.*f(x); 
% Plot nullrule
plot(e);
xlabel('Graad');
ylabel('Error');

% Bepaal de reductiefactoren
for j = 1: length(e)-1
    if (e(j+1)>10^(-5))
        r(j) = e(j)/e(j+1);
    end    
end
rm = max(r);

% Vermijdt het fase effect
for j = 1: length(e)/2
    E(j)=sqrt(e(2j-1)^2+e(2j)^2);
end

% Bepaal de reductiefactoren
for j = 1: length(E)-1
    if (E(j+1)>10^(-5))
        R(j) = E(j)/E(j+1);
    end    
end

Rm = max(R);
figure();
plot(r);
hold on;
plot(R);
xlabel('Graad');
ylabel('Reductiefactoren');
legend('r', 'R');
title('Verband tussen reductiefactoren en fase-effecten');
end

