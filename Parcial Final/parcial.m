filename = 'vidas.txt';
[Vida,~]=importdata(filename);
%% 

%Graficar distribuciones empiricas
clf

figure(1)
CM = parula(17);
for j = 1:13
    
    [F,t] = ecdf(Vida(:,j));
    plot(t,F, 'color', CM(j,:))
   
    hold on
    xlabel('t')
    ylabel('ECDF')
    legend

end 


% si queremos que el tiempo sea mayor en promedio, buscamos la curva que
% este mas abajo ya que el area que encontramos es la que esta arriba de la
% curva.

%los ultimos anos tienen mayor tiempo de vida en promedio

%%
%Bandas de confianza para la empirica

Confianza = 0.9;
alpha = 1 - Confianza;
epsilon = sqrt(1/(2*length(Vida))*log(2/alpha));

%minimo - bandas minimo 

[F, t] = ecdf(Vida(:,1));
upper_band_min = zeros(length(F),1);lower_band_min=zeros(length(F),1);
for j = 1:length(F)
    upper_band_min(j) = min(F(j)+epsilon,1);
    lower_band_min(j) = max(F(j)-epsilon,0);
end


clf 
figure(8)
%plotting coldest
plot(t,F, 'black')
hold on 
plot(t,upper_band_min,'blue')
hold on 
plot(t, lower_band_min,'blue')

legend('Empirica','Banda Superior','Banda Inferior')
xlabel('t')
ylabel('ECDF')
title('Bandas de confianza para el primer ano')
hold off

%%

[F, t] = ecdf(Vida(:,end));
mu = 13;
p = expcdf(t,mu);
plot(t,F, 'black')
hold on
plot(t,p, 'blue')

legend('Empirica ultimo ano','Exponencial media 13')
xlabel('t')
ylabel('CDF')
title('Ultimo ano comparado con exponencial de media 13')
hold off


%%
%Densidades estimadas y aproximadas
[f,xi] = ksdensity(Vida(:,end));
[f1,xj] = ksdensity(Vida(:,end),'Kernel','epanechnikov');
[f2,xk] = ksdensity(Vida(:,end),'Kernel',@exppdf);
p = exppdf(t,mu);
plot(xi,f)
hold on
plot(xj,f1)
plot(xk,f2)
plot(t,p, 'LineWidth',1)
legend('Kernel Gaussiano','Kernel Epanechnikov','Kernel Exponencial','Exponencial media 13')
title('Comparacion de las densidades estimada y aproximadas')
xlabel('t')
ylabel('PDF')


%%
%Glivenko Cantelli
diff = zeros(length(Vida(:,end))-1,1);
mu = 13;
for i = 2:length(Vida(:,end))
    [F, t] = ecdf(Vida(1:i,end));
    p = expcdf(t,mu);
    diff(i-1) = max(abs(F-p));
end

plot(diff,'color','#A2142F')
title("Glivenko Cantelli theorem demonstration")
xlabel("Sample size (n)")
ylabel("Maximum Error")

%%
%Intervalo de confianza  bootstrap

boot = bootstrp(10000,@max,Vida(:,1));
max(Vida(:,1))
%confidence interval
CIB = [prctile(boot,2.5) prctile(boot,97.5)] 
%clf
%hist(boot);

jack = jackknife(@max,Vida(:,1));
n = length(Vida(:,1));
jbias = (n-1)* (mean(jack)-max(Vida(:,1))) % jackknife bias estimate


var(jack)

%%
filename = 'costos.txt';
[Costos,~]=importdata(filename);

%%
%Regresion
t = Vida(:,1);

linmod = fitlm(t,Costos);

NadarayaWatson = ksr(t,Costos)
plot(NadarayaWatson.x, NadarayaWatson.f,'black')
hold on 
plot(t,Costos,'.b')
legend({'Regresion de Nadaraya Watson','Datos'})

xlabel('Primer ano')
ylabel('Costos')

title('Regresion de Nadaraya Watson')

Nada = interp1(NadarayaWatson.x,NadarayaWatson.f,t);

OrdinaryR2 = 1- sum((Costos-Nada).^2)/sum((Costos-mean(Costos)).^2);
AjustedR2 =  OrdinaryR2*((length(Costos)-1)/(length(Costos)-2));


%%
%Test de rangos
%Pasamos el test de MannWhitney
N = size(Vida,2); n = size(Vida,1);
MannWhitneyp = zeros(N,N); MannWhitneyh = zeros(N,N);
for i = 1:N
    Xi = Vida(:,i);
    for j = 1:N
        Xj = Vida(:,j);
        [p,h] = ranksum(Xi,Xj);
        MannWhitneyp(i,j) = p; MannWhitneyh(i,j) = h;
    end
end
I = MannWhitneyp <0.05
%Observamos que ninguno no pasa el test, entonces vamos a acortar la
%ventana de tiempo.
heatmap(double(I))
colormap spring

%%
%Test de homogeneidad
X  = Vida(1:400,:); Y  = Vida(end-399:end,:);
ddplot2(X,Y)
%%
[LM, PD1, PD2] = TestDepthTukey(X, Y);
