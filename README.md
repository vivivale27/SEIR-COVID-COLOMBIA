function colombia_covid_model

clear all;
close all;

global Tinf Ntot Tinc Tsev Tdeath K1 K2 K3 K4 pa pm ps pf mix_matrix pop_dist R0 num_variants

% Períodos característicos
Tinc = 4.9;  % periodo de incubación
Tinf = 5;    % periodo infeccioso
Tsev = 22.6; % duración de casos severos (hospitalización)
Tdeath = 9;  % duración hasta fallecimiento tras hospitalización

% Números de etapas en distribución Gamma
K1 = 2;  % etapas de incubación
K2 = 4;  % etapas infecciosas
K3 = 24; % etapas de hospitalización
K4 = 2;  % etapas hasta fallecimiento

% Distribución poblacional de Colombia (16 grupos etarios de 5 años)
% Datos del DANE 2020 (en miles de personas)
pop_dist = [4342 4264 4303 4345 4400 4325 4088 3860 3460 2939 2333 1954 1530 1080 687 864];
pop_dist = pop_dist * 1000; % convertir a número de personas
Ntot = sum(pop_dist);

% Cargar o definir manualmente parámetros epidemiológicos
% En un caso real estos se cargarían desde archivos de datos
psymp = ones(16,1) * 0.7;      % proporción sintomática
pasymp = ones(16,1) * 0.3;     % proporción asintomática
hosp_rates = ones(16,1) * 0.15; % tasas de hospitalización
mort_rates = ones(16,1) * 0.02; % tasas de mortalidad

% Ajustar tasas por edad
for i = 1:16
    factor_edad = 1 + (i-1)*0.1; % Factor que aumenta con la edad
    if i > 12
        hosp_rates(i) = min(0.4, hosp_rates(i) * factor_edad);
        mort_rates(i) = min(0.15, mort_rates(i) * factor_edad);
    end
end

pa = pasymp; % proporción asintomática
ph = hosp_rates; % proporción que requiere hospitalización
pm = psymp - ph; % proporción sintomática que no requiere hospitalización
pf = mort_rates; % proporción que fallece
ps = ph - pf; % proporción que requiere hospitalización pero se recupera

% Generar matrices de contacto simplificadas para Colombia
% Normalmente se cargarían desde archivos
school = generate_contact_matrix(0.3, 16, 'school');
work = generate_contact_matrix(0.5, 16, 'work');
other = generate_contact_matrix(0.2, 16, 'other');
home = generate_contact_matrix(0.4, 16, 'home');

% Combinar matrices de contacto
a1 = 1; % escuela
a2 = 1; % trabajo
a3 = 1; % otros
a4 = 1; % hogar
mix_matrix = a1*school + a2*work + a3*other + a4*home;

% Número reproductivo básico
R0 = 2.5;

% Simulación del modelo SEIR
IC0 = zeros(1, (K1+K2+K3+K4+5)*16);
IC0(1:16) = pop_dist; % Inicialmente toda la población es susceptible
IC0(8) = IC0(8) - 100; % Introducir 100 casos en el grupo etario 8 (30-35 años)
IC0(16*1 + 8) = 100; % Estos casos pasan a ser expuestos

% Configurar opciones para el solver
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Resolver el sistema de ecuaciones diferenciales
[t, x] = ode45(@seir_age, [0 100], IC0, opts);

% Procesar resultados para visualización
total_cases = zeros(length(t), 1);
total_deaths = zeros(length(t), 1);
for i = 1:length(t)
    % Sumar todos los que han pasado por la infección (I, Ra, Rm, Rs, IF, D)
    total_cases(i) = sum(x(i, 16*(K1+1)+1:end));
    % Sumar fallecidos
    total_deaths(i) = sum(x(i, 16*(K1+K2+K3+K4+4)+1:16*(K1+K2+K3+K4+5)));
end

% Visualizar resultados
figure(1);
subplot(2,1,1);
plot(t, total_cases, 'b-', 'LineWidth', 2);
xlabel('Días');
ylabel('Casos totales');
title('Proyección de casos COVID-19 en Colombia');
grid on;

subplot(2,1,2);
plot(t, total_deaths, 'r-', 'LineWidth', 2);
xlabel('Días');
ylabel('Muertes totales');
title('Proyección de muertes COVID-19 en Colombia');
grid on;

disp('Simulación completada.');

end

% Función auxiliar para generar matrices de contacto simplificadas
function matrix = generate_contact_matrix(base_value, num_groups, type)
    matrix = ones(num_groups) * base_value;
    
    if strcmp(type, 'school')
        % Mayor contacto entre grupos jóvenes
        for i = 1:5
            for j = 1:5
                matrix(i,j) = matrix(i,j) * 3;
            end
        end
    elseif strcmp(type, 'work')
        % Mayor contacto entre adultos (grupos 4-12)
        for i = 4:12
            for j = 4:12
                matrix(i,j) = matrix(i,j) * 2;
            end
        end
    elseif strcmp(type, 'home')
        % Contacto intergeneracional
        for i = 1:num_groups
            for j = 1:num_groups
                if abs(i-j) <= 2 || (i <= 5 && j >= 13) || (j <= 5 && i >= 13)
                    matrix(i,j) = matrix(i,j) * 1.5;
                end
            end
        end
    end
end

% Implementación del modelo SEIR con estructura por edad
function dxdt = seir_age(t, x)
    global Tinf Tinc Tsev Tdeath K1 K2 K3 K4 R0 pa pf ps pm mix_matrix

    % Número de grupos etarios
    num_age_groups = 16;
    
    % Inicializar el vector de derivadas
    dxdt = zeros(size(x));
    
    % Calcular población total por grupo etario
    NN = zeros(num_age_groups, 1);
    for i = 1:num_age_groups
        for j = 0:(K1+K2+K3+K4+4)
            idx = j*num_age_groups + i;
            if idx <= length(x)
                NN(i) = NN(i) + x(idx);
            end
        end
    end
    
    % Calcular valor propio dominante de la matriz de contacto
    r1 = max(real(eig(mix_matrix)));
    
    % Calcular infectados actuales por grupo etario
    Icur = zeros(num_age_groups, 1);
    for i = 1:num_age_groups
        for j = 1:K2
            idx = (K1+j)*num_age_groups + i;
            if idx <= length(x)
                Icur(i) = Icur(i) + x(idx);
            end
        end
    end
    
    % Ecuaciones para susceptibles (S)
    for i = 1:num_age_groups
        idx_S = i;
        sum_force = 0;
        
        for j = 1:num_age_groups
            sum_force = sum_force + mix_matrix(i,j) * Icur(j) / (Tinf * NN(j));
        end
        
        dxdt(idx_S) = -R0 * x(idx_S) * sum_force / r1;
    end
    
    % Ecuaciones para expuestos (E1)
    for i = 1:num_age_groups
        idx_S = i;
        idx_E1 = num_age_groups + i;
        
        dxdt(idx_E1) = -dxdt(idx_S) - K1 * x(idx_E1) / Tinc;
    end
    
    % Ecuaciones para el resto de compartimentos expuestos (E2 a E_K1)
    for i = 1:num_age_groups
        for j = 1:(K1-1)
            idx_Ej = (j)*num_age_groups + i;
            idx_Ej_next = (j+1)*num_age_groups + i;
            
            dxdt(idx_Ej_next) = K1 * x(idx_Ej) / Tinc - K1 * x(idx_Ej_next) / Tinc;
        end
    end
    
    % Ecuaciones para infectados (I1)
    for i = 1:num_age_groups
        idx_EK1 = K1*num_age_groups + i;
        idx_I1 = (K1+1)*num_age_groups + i;
        
        dxdt(idx_I1) = K1 * x(idx_EK1) / Tinc - K2 * x(idx_I1) / Tinf;
    end
    
    % Ecuaciones para el resto de compartimentos infectados (I2 a I_K2)
    for i = 1:num_age_groups
        for j = 1:(K2-1)
            idx_Ij = (K1+j)*num_age_groups + i;
            idx_Ij_next = (K1+j+1)*num_age_groups + i;
            
            dxdt(idx_Ij_next) = K2 * x(idx_Ij) / Tinf - K2 * x(idx_Ij_next) / Tinf;
        end
    end
    
    % Ecuaciones para recuperados asintomáticos (Ra)
    for i = 1:num_age_groups
        idx_IK2 = (K1+K2)*num_age_groups + i;
        idx_Ra = (K1+K2+1)*num_age_groups + i;
        
        dxdt(idx_Ra) = K2 * pa(i) * x(idx_IK2) / Tinf;
    end
    
    % Ecuaciones para recuperados con síntomas leves (Rm)
    for i = 1:num_age_groups
        idx_IK2 = (K1+K2)*num_age_groups + i;
        idx_Rm = (K1+K2+2)*num_age_groups + i;
        
        dxdt(idx_Rm) = K2 * pm(i) * x(idx_IK2) / Tinf;
    end
    
    % Ecuaciones para infectados severos (IS1)
    for i = 1:num_age_groups
        idx_IK2 = (K1+K2)*num_age_groups + i;
        idx_IS1 = (K1+K2+3)*num_age_groups + i;
        
        dxdt(idx_IS1) = K2 * ps(i) * x(idx_IK2) / Tinf - K3 * x(idx_IS1) / Tsev;
    end
    
    % Ecuaciones para el resto de compartimentos de infectados severos (IS2 a IS_K3)
    for i = 1:num_age_groups
        for j = 1:(K3-1)
            idx_ISj = (K1+K2+2+j)*num_age_groups + i;
            idx_ISj_next = (K1+K2+3+j)*num_age_groups + i;
            
            if idx_ISj <= length(x) && idx_ISj_next <= length(x)
                dxdt(idx_ISj_next) = K3 * x(idx_ISj) / Tsev - K3 * x(idx_ISj_next) / Tsev;
            end
        end
    end
    
    % Ecuaciones para recuperados severos (Rs)
    for i = 1:num_age_groups
        idx_ISK3 = (K1+K2+K3+2)*num_age_groups + i;
        idx_Rs = (K1+K2+K3+3)*num_age_groups + i;
        
        if idx_ISK3 <= length(x) && idx_Rs <= length(x)
            dxdt(idx_Rs) = K3 * x(idx_ISK3) / Tsev;
        end
    end
    
    % Ecuaciones para infectados fatales (IF1)
    for i = 1:num_age_groups
        idx_IK2 = (K1+K2)*num_age_groups + i;
        idx_IF1 = (K1+K2+K3+4)*num_age_groups + i;
        
        if idx_IK2 <= length(x) && idx_IF1 <= length(x)
            dxdt(idx_IF1) = K2 * pf(i) * x(idx_IK2) / Tinf - K4 * x(idx_IF1) / Tdeath;
        end
    end
    
    % Ecuaciones para el resto de compartimentos de infectados fatales (IF2 a IF_K4)
    for i = 1:num_age_groups
        for j = 1:(K4-1)
            idx_IFj = (K1+K2+K3+3+j)*num_age_groups + i;
            idx_IFj_next = (K1+K2+K3+4+j)*num_age_groups + i;
            
            if idx_IFj <= length(x) && idx_IFj_next <= length(x)
                dxdt(idx_IFj_next) = K4 * x(idx_IFj) / Tdeath - K4 * x(idx_IFj_next) / Tdeath;
            end
        end
    end
    
    % Ecuaciones para fallecidos (D)
    for i = 1:num_age_groups
        idx_IFK4 = (K1+K2+K3+K4+3)*num_age_groups + i;
        idx_D = (K1+K2+K3+K4+4)*num_age_groups + i;
        
        if idx_IFK4 <= length(x) && idx_D <= length(x)
            dxdt(idx_D) = K4 * x(idx_IFK4) / Tdeath;
        end
    end
end
