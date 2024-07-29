function Foods=generate_random_solution
% global a1_min a1_max a2_min a2_max a3_min a3_max
global a_min a_max order
% Foods(1)=rand(1)*(a1_max-a1_min)+a1_min;
% Foods(2)=rand(1)*(a2_max-a2_min)+a2_min;
% Foods(3)=rand(1)*(a3_max-a3_min)+a3_min;
for i=1:order+1
    Foods(i)=rand(1)*(a_max-a_min)+a_min;
end
end