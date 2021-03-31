clear; clc

error_rate = 0.01;
%% Assuming Symetric in One Axis
material_stress_yield = 350;
material_strain_yield = 0.0018;
material_elasticity_modulus = 200000; %N/mm^2
% Section Data
section_width = [200; 11; 200];
section_width_coordinates = [0 17; 17 583; 583 600];

% n-pieces division
N = 200;

% Load Data
Axial_Load = 0;
N_geometrics = [];
%% Generate N-pieces Geometric Data and Stress-Strain Matrix
section_height = max(max(section_width_coordinates));
element_height = section_height / N;
centroid = element_height / 2; %% From the element top/bottom.
% N_geometrics = [%Y_start, Y_end, centroid_relative, centroid_abs,
% element_width]
N_geometrics = [0 0+element_height, centroid, centroid, section_width(find(section_width_coordinates(:,1)<=centroid & section_width_coordinates(:,2) > centroid,1,'first'))];
for n=2:N
    N_geometrics = [N_geometrics; N_geometrics(n-1,2), N_geometrics(n-1,2)+element_height, centroid, N_geometrics(n-1,2)+centroid, section_width(find(section_width_coordinates(:,1)<=N_geometrics(n-1,2)+centroid & section_width_coordinates(:,2) > N_geometrics(n-1,2)+centroid,1,'first'))];
end
N_stressstrain = zeros(N, 2);


%% Start Lamina Analysis
% assume Neutral Axis Location
curvature_array = [0:1:1]*1000;
NA_height = section_height / 2;
%% NA_height generator
NA_heights = 1:section_height; %% A more exotic algorithm is required if the NA is not integer.

NA_height_iterator = 1;
analysis_correct = false;
iteration_data = [];
curvature_array = [0:0.0000001:0.00005];
curvature_na_data = [];
%% Calculate Moment at TOTAL yield Condition
plastic_moment = 0;
plastic_curvature = 0;
yield_moment = 0;
yield_curvature = 0;
while true
    if NA_height_iterator <= length(NA_heights)
        NA_height = NA_heights(NA_height_iterator);
    end
    for curvature_iteration=1:length(curvature_array)   
        analysis_correct = false;
        curvature = curvature_array(curvature_iteration);
        N_geometrics(:,6) = [N_geometrics(:,4)-NA_height];
        N_stressstrain(:,2) = (N_geometrics(:,6))*curvature;
        N_stressstrain(:,1) = N_stressstrain(:,2) *material_elasticity_modulus;
        N_stressstrain(N_stressstrain(:,2)>=material_strain_yield,1) = material_stress_yield;
        N_stressstrain(N_stressstrain(:,2)<=-material_strain_yield,1) = -material_stress_yield;
        section_load = sum(N_stressstrain(:,1).*element_height.*N_geometrics(:,5));
        section_moment = sum(N_stressstrain(:,1).*element_height.*N_geometrics(:,5).*(N_geometrics(:,6)));
        if (abs(section_load - Axial_Load) < error_rate)
                analysis_correct = true;
                if (((abs(N_stressstrain(1,2)) >= material_strain_yield) |  (abs(N_stressstrain(end,2)) >= material_strain_yield)) & yield_moment ==0)
                      yield_moment = section_moment;
                     yield_curvature = curvature;
                
                end
                curvature_na_data = [curvature_na_data; NA_height, curvature, section_load, section_moment, 1];
                continue 
        end
        iteration_data = [iteration_data; NA_height, curvature,section_load, section_moment, 0];
        NA_height_iterator = NA_height_iterator+1;
        break
        
    end
    if (analysis_correct == true)
        str = "NA is "+ NA_height
        break
    end

end

curvature_na_data(:,7) = curvature_na_data(:,4)/yield_moment;
curvature_na_data(:,6) = curvature_na_data(:,2)/yield_curvature;
plot(curvature_na_data(:,6),curvature_na_data(:,7))
xlim([0 4])
ylim([0 1.2])
grid on
