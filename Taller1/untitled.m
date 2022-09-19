%% Punto 1
%Graficar empíricas
x = 1:35;
names = string(x);
[~,n] = size(temp);
for i = 1:n
    cdfplot(temp(:,i));
    hold on
end
hold off
legend(names)

%Estimar media y media empírica
m = zeros(35,35);
ms = zeros(35,35);
media = mean(temp); %Media de cada año

%Tabla doble entrada media normal
for i=1:35
    for j =1:35
        if media(i)<media(j)
            m(i,j) = 0;
        end
        if media(i)==media(j)
            m(i,j) = 0.5;
        end
        if media(i)>media(j)
            m(i,j) = 1;
        end
    end
end
heatmap(m)

%Tabla doble entrada orden estocástico

empiricas = zeros(1,35);
for i=1:35
    empiricas(i) = ecdf(temp(:,i));
end



%% Punto 2
media_est = zeros(35,1);

%Estimar media empírica
for i=1:35
    temps = temp(:,i);
    [F,t] = ecdf(temps);
    idx1 = t(find(t>=0));
    idx2 = t(find(t<0));
    Fq1 = F(find(t>=0));
    Fq2 = F(find(t<0));
    estimator = trapz(idx1,1-Fq1) - trapz(idx2,Fq2);
    media_est(i) = estimator;
end
clf
plot(1:35, media_est, 'o')
hold on
plot(1:35, media, 'o')
title("Estimator plug-in vs Media")
xlabel("Año")
ylabel("Temperatura media estimada")
legend('Estimador','Media')

%% Punto 3

%Encontrar caluroso y frío
media = means;
maxt = 1;
mint = 1;
for i=1:35
    if media(i)>media(maxt)
        maxt = i;
    end
    if media(i)<media(mint)
        mint = i;
    end
end
frio = temp(:,mint);
calor = temp(:, maxt);
%Graficar empíricas con bandas de confianza 
ecdf(calor,'Bounds','on')
hold on
ecdf(frio,'Bounds','on')
hold on
legend('Max temp','','', 'Min temp')
%% Punto 12
%Intervalo de confianza para la temperatura máxima
boot = bootstrp(1000,@max,minYear);
sub = prctile(boot,2.5);
sup = prctile(boot,97.5);
bootci(1000,@max,minYear)
ci = [sub sup]
%Sesgo del parámetro por jackknife
n = size(minYear);
jack = jackknife(@max,minYear);
jbias = (n-1)*(mean(jack)-max(minYear))
%Varianza
varMax = var(boot)

%% Punto 15
clf
frio = minYear
calor = maxYear
dist = mahal(calor, frio);
p = prctile(dist,95);
idx = find(dist>p);
plot(calor, frio, 'o')
hold on
plot(Temps(idx,maxt), Temps(idx,mint), 'o')
title("Detección de outliers Mahalanobis clásica")
xlabel("Temperaturas más calurosas en media")
ylabel("Temperatura más frías en media")
legend('','Outliers')
%%
clf
x = [calor frio];
cova = cov(x);
x_minus_mu = x-media(:,[maxt mint]);
[P,R,C] = equilibrate(cova);
cova = R*P*cova*C;
invcov = inv(cova);
left_term = x_minus_mu*invcov;
mahal_nueva = left_term* x_minus_mu';
dist = diag(mahal_nueva);
p = prctile(dist,95);
idx = find(dist>p);
plot(calor, frio, 'o')
hold on
plot(Temps(idx,maxt), Temps(idx,mint), 'o')
