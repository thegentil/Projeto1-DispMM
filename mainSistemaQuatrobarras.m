function [] = mainSistemaQuatrobarras()

% Definindo os inputs

theta0 = 30*pi/180;
vin = 0.25; % m/s
l = 0.3; % m
g = 9.81; % m/s^2
m = 1; % Massa da barra 1
m2 = 0.25; % Massa das barras 2,3,4 e 5
mr = 0.1; % Massa do rolete
Gy = 300; %Força
lvar = linspace(0,1.5,100);
m2var = linspace(0.5,10,100);

tmax = l*cos(theta0)/vin;

% Criando a série de tempo (linspace)

tempo = linspace(0,tmax,120);

% Equacionando theta

x = cos(theta0)-vin.*tempo/l

theta = acos(x);
theta1 = (90*pi)/180 - theta;

% CINEMÁTICA

% Equacionando as velocidades 

Vbx = vin;
Vby = vin.*cot(theta);

Vg = 2*vin.*cot(theta); % Vgy

wbf = (-vin)./(l.*sin(theta));
wgb = (vin)./(l.*sin(theta));

% Equacionando as acelerações 

alphabf = ((vin./l)^2).*(cos(theta)./(sin(theta).^3));
alphagb = -(wgb.^2).*cot(theta);


Ab = (-alphabf).*l.*cos(theta) - (wbf.^2).*l.*sin(theta); % Aby

Ag = Ab + alphagb.*l.*cos(theta) - (wgb.^2).*l.*sin(theta); 

Acm4x = (-0.5).*alphagb.*l.*sin(theta) - 0.5.*(wgb.^2).*l.*cos(theta);

Acm4y = Ab + 0.5.*alphagb.*l.*cos(theta) - 0.5.*(wgb.^2).*l.*sin(theta);

Acm5x = (-0.5).*alphabf.*l.*sin(theta) + 0.5.*(wbf.^2).*l.*cos(theta);

Acm5y = (-0.5).*alphabf.*l.*cos(theta) - 0.5.*(wbf.^2).*l.*sin(theta);

Ary = Ab;

Arx = 0;

% DINÂMICA VARIANDO O PARÂMETRO l (comprimento das barras 4 e 5) considerando a
% angulação da barra de 30 graus constante (theta0) e a massa da
% barra perfeita.

wbflvar = (-vin)./(lvar.*sin(theta0));
alphabflvar = ((vin./lvar).^2).*(cos(theta0)./(sin(theta0).^3));
wgblvar = (vin)./(lvar.*sin(theta0));
alphagblvar = -(wgblvar.^2).*cot(theta0);
Ablvar = (-alphabflvar).*lvar.*cos(theta0) - (wbflvar.^2).*lvar.*sin(theta0);
Acm4xlvar = (-0.5).*alphagblvar.*lvar.*sin(theta0) - 0.5.*(wgblvar.^2).*lvar.*cos(theta0);
Acm4ylvar = Ablvar + 0.5.*alphagblvar.*lvar.*cos(theta0) - 0.5.*(wgblvar.^2).*lvar.*sin(theta0);
Arylvar = Ablvar;
Acm5xlvar = (-0.5).*alphabflvar.*lvar.*sin(theta0) + 0.5.*(wbflvar.^2).*lvar.*cos(theta0);
Acm5ylvar = (-0.5).*alphabflvar.*lvar.*cos(theta0) - 0.5.*(wbflvar.^2).*lvar.*sin(theta0);
Arxlvar = 0;
Gylvar = 300;

Fy4lvar = m2.*(Acm4ylvar + g) + Gylvar;
Fy5lvar = mr.*(Arylvar + g) + Fy4lvar;
Fx5lvar = ((m2.*lvar.*alphagblvar)/(3.*sin(theta0))) - ((Fy5lvar.*cos(theta0))/sin(theta0));
Oxlvar = m2.*Acm5xlvar + Fx5lvar;
Oylvar = m2.*(Acm5ylvar + g) + Fy5lvar;
Gxlvar = ((- m2.*lvar.*alphabflvar)./(12.*sin(theta0))) - ((Fy4lvar.*cos(theta0))./(2.*sin(theta0))) - ((Gylvar.*cos(theta0))./(2.*sin(theta0))) + m2.*Acm4xlvar;
Fx4lvar = Gxlvar - m2.*Acm4xlvar;
Fmlvar = -(Fx5lvar + Fx4lvar + mr.*Arxlvar).*2;

% DINÂMICA VARIANDO O PARÂMETRO m (massa das barras 4 e 5) considerando a
% angulação da barra de 30 graus constante (theta0) e o comprimento de
% barra perfeito.
m2

wbfmvar = (-vin)./(l.*sin(theta0));
alphabfmvar = ((vin./l).^2).*(cos(theta0)./(sin(theta0).^3));
wgbmvar = (vin)./(l.*sin(theta0));
alphagbmvar = -(wgbmvar.^2).*cot(theta0);
Abmvar = (-alphabfmvar).*l.*cos(theta0) - (wbfmvar.^2).*l.*sin(theta0);
Acm4xmvar = (-0.5).*alphagbmvar.*l.*sin(theta0) - 0.5.*(wgbmvar.^2).*l.*cos(theta0);
Acm4ymvar = Abmvar + 0.5.*alphagbmvar.*l.*cos(theta0) - 0.5.*(wgbmvar.^2).*l.*sin(theta0);
Arymvar = Abmvar;
Acm5xmvar = (-0.5).*alphabfmvar.*l.*sin(theta0) + 0.5.*(wbfmvar.^2).*l.*cos(theta0);
Acm5ymvar = (-0.5).*alphabfmvar.*l.*cos(theta0) - 0.5.*(wbfmvar.^2).*l.*sin(theta0);
Arxmvar = 0;
Gymvar = 300;

Fy4mvar = m2var.*(Acm4ymvar + g) + Gymvar;
Fy5mvar = mr.*(Arymvar + g) + Fy4mvar;
Fx5mvar = ((m2var.*l.*alphagbmvar)/(3.*sin(theta0))) - ((Fy5mvar.*cos(theta0))/sin(theta0));
Oxmvar = m2var.*Acm5xmvar + Fx5mvar;
Oymvar = m2var.*(Acm5ymvar + g) + Fy5mvar;
Gxmvar = ((- m2var.*l.*alphabfmvar)./(12.*sin(theta0))) - ((Fy4mvar.*cos(theta0))./(2.*sin(theta0))) - ((Gymvar.*cos(theta0))./(2.*sin(theta0))) + m2.*Acm4xmvar;
Fx4mvar = Gxmvar - m2var.*Acm4xmvar;
Fmmvar = -(Fx5mvar + Fx4mvar + mr.*Arxmvar).*2;

% DINÂMICA COM OS PARÂMETROS PERFEITOS

Fy4 = m2.*(Acm4y + g) + Gy;
Fy5 = mr.*(Ary + g) + Fy4;
Fx5 = ((m2.*l.*alphagb)/(3.*sin(theta))) - ((Fy5.*cos(theta))/sin(theta));
Ox = m2.*Acm5x + Fx5;
Oy = m2.*(Acm5y + g) + Fy5;
Gx = ((- m2.*l.*alphabf)/(12*sin(theta))) - ((Fy4.*cos(theta))/(2.*sin(theta))) - ((Gy.*cos(theta))/(2.*sin(theta))) + m2.*Acm4x;
Fx4 = Gx - m2.*Acm4x;
Fm = -(Fx5 + Fx4 + mr.*Arx).*2;

F4 = (Fy4.^2 + Fx4.^2).^(1/2);
F5 = (Fx5.^2 + Fy5.^2).^(1/2);
O = (Ox.^2 + Oy.^2).^(1/2);
G = (Gx.^2 + Gy.^2).^(1/2);
Gplot = 2.*G;


% Janela 1: wab x t | wca x t

figure(1)

whitebg('k')

subplot (1,2,1)
plot(tempo,wbf,'r','LineWidth',2)
title('Velocidade Angular AB x Tempo')
xlabel('Tempo (s)')
ylabel('Wab (rad/s)')
grid
axis tight

subplot(1,2,2)
plot(tempo,wgb,'w','LineWidth',2)
title('Velocidade Angular CA x Tempo')
xlabel('Tempo (s)')
ylabel('Wca (rad/s)')
grid
axis tight    

% Janela 2: Vay x t | Vay x Theta

figure(2)

whitebg('k')

subplot (1,2,1)
plot(tempo,Vby,'r','LineWidth',2)
title('Velocidade ponto A x Tempo')
xlabel('Tempo (s)')
ylabel('Vay (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Vby,'w','LineWidth',2)
title('Velocidade ponto A x Ângulo')
xlabel('Theta (graus)')
ylabel('Vay (m/s)')
grid
axis tight

% Janela 3: wbf x t | wbf x Theta

figure(3)

whitebg('k')

subplot (1,2,1)
plot(tempo,wbf, 'r','LineWidth',2)
title('Veloc. Âng. barra AB x Tempo')
xlabel('Tempo (s)')
ylabel('W barra AB (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,wbf,'w','LineWidth',2)
title('Veloc. Âng. barra AB x Ângulo')
xlabel('Theta (graus)')
ylabel('W barra AB (m/s)')
grid
axis tight

% Janela 4: wgb x t | wgb x Theta

figure(4)

whitebg('k')

subplot (1,2,1)
plot(tempo,wgb,'r','LineWidth',2)
title('Veloc. Âng. barra AC x Tempo')
xlabel('Tempo (s)')
ylabel('W barra AC (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,wgb, 'w','LineWidth',2)
title('Veloc. Âng. barra AC x Ângulo')
xlabel('Theta (graus)')
ylabel('W barra AC (m/s)')
grid
axis tight

% Janela 5: Vcy x t | Vcy x Theta

figure(5)

whitebg('k')

subplot (1,2,1)
plot(tempo,Vg,'r','LineWidth',2)
title('Velocidade ponto C x Tempo')
xlabel('Tempo (s)')
ylabel('Vcy (m/s)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Vg, 'w','LineWidth',2)
title('Velocidade ponto C x Ângulo')
xlabel('Theta (graus)')
ylabel('Vcy (m/s)')
grid
axis tight


% Janela 6: Ag x t | Ag x Theta

figure (6)

whitebg('k')

subplot (1,2,1)
plot(tempo,Ag,'r','LineWidth',2)
title('Aceleração ponto C x Tempo')
xlabel('Tempo (s)')
ylabel('Acy (m/s^2)')
grid
axis tight

subplot(1,2,2)
plot(theta*180/pi,Ag,'w','LineWidth',2)
title('Aceleração ponto C x Ângulo')
xlabel('Theta (graus)')
ylabel('Acy (m/s^2)')
grid
axis tight

% Janela 7: Ab x t | Ab x Theta

figure (7)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Ab,'r','LineWidth',2)
title('Aceleração ponto A x Tempo')
xlabel('Tempo (s)')
ylabel('Ay (m/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,Ab,'w','LineWidth',2)
title('Aceleração ponto A x Ângulo')
xlabel('Theta (graus)')
ylabel('Ay (m/s^2)')
grid
axis tight

figure (8)

whitebg('k') 

subplot (1,2,1)
plot(tempo,alphabf,'r','LineWidth',2)
title('Alpha barra AB x Tempo')
xlabel('Tempo (s)')
ylabel('Alpha barra AB (rad/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,alphabf,'w','LineWidth',2)
title('Alpha barra AB x Ângulo')
xlabel('Theta (graus)')
ylabel('Alpha barra AB (rad/s^2)')
grid
axis tight

figure (9)

whitebg('k') 

subplot (1,2,1)
plot(tempo,alphagb,'r','LineWidth',2)
title('Alpha barra AC x Tempo')
xlabel('Tempo (s)')
ylabel('Alpha barra AC (rad/s^2)')
grid
axis tight

subplot (1,2,2)
plot(theta*180/pi,alphagb,'w','LineWidth',2)
title('Aceleração barra AC x Ângulo')
xlabel('Theta (graus)')
ylabel('Alpha barra AC (rad/s^2)')
grid
axis tight

figure (10)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Fy4,'r','LineWidth',2)
title('F4y x Tempo')
xlabel('Tempo')
ylabel('F4y (N)')
grid
axis tight

subplot (1,2,2)
plot(tempo,Fx4,'w','LineWidth',2)
title('F4x x Tempo')
xlabel('Tempo')
ylabel('F4x (N)')
grid
axis tight

figure (11)

whitebg('k') 

plot(tempo,Fy5,'r','LineWidth',2)
title('F5y x Tempo')
xlabel('Tempo')
ylabel('F5y (N)')
grid
axis tight


figure (12)

whitebg('k') 

plot(tempo,Fm,'r','LineWidth',2)
title('Força de entrada x Tempo')
xlabel('Tempo (s)')
ylabel('Força de entrada (N)')
grid
axis tight


figure (13)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Gplot,'r','LineWidth',2)
title('Força de saída x Tempo')
xlabel('Tempo (s)')
ylabel('Força de saída (N)')
grid
axis tight

subplot (1,2,2)
plot(tempo,Fm,'w','LineWidth',2)
title('Força de entrada x Tempo')
xlabel('Tempo (s)')
ylabel('Força de entrada (N)')
grid
axis tight

figure (14)

whitebg('k') 

plot(lvar,Fmlvar,'r','LineWidth',2)
title('Força de entrada x Comprimento das barras 4 e 5')
xlabel('Comprimento da barra (m)')
ylabel('Força de entrada (N)')
grid
axis tight

figure (15)

whitebg('k') 

plot(m2var,Fmmvar,'r','LineWidth',2)
title('Força de entrada x Massa das barras 4 e 5')
xlabel('Massa das barras (kg)')
ylabel('Força de entrada (N)')
grid
axis tight

figure (16)

whitebg('k') 

subplot (1,2,1)
plot(tempo,alphabf,'r','LineWidth',2)
title('alpha barra AB x Tempo')
xlabel('Tempo (s)')
ylabel('alpha barra AB (rad/s^2)')
grid
axis tight

subplot (1,2,2)
plot(tempo,alphagb,'w','LineWidth',2)
title('alpha barra CA x Tempo')
xlabel('Tempo (s)')
ylabel('alpha barra CA (rad/s^2)')
grid
axis tight

figure (17)

whitebg('k') 

subplot (1,2,1)
plot(tempo,Oy,'r','LineWidth',2)
title('FBy x Tempo')
xlabel('Tempo')
ylabel('FBy (N)')
grid
axis tight

subplot (1,2,2)
plot(tempo,Ox,'w','LineWidth',2)
title('FBx x Tempo')
xlabel('Tempo')
ylabel('FBx (N)')
grid
axis tight

for i = 1:119
    integral(i) = Ag(i).*(tempo(i + 1) - tempo(i));
end

aef = cumsum(integral)
tempo
tempo*vin

end