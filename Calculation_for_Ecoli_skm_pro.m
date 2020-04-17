% initialization
    % cobra_open
    %initCobraToolbox
    
% Read yeast model
%Inputfile = './Dataset/Input/input_SA5_pTH-aroG_fbr-ppsA-tktA.tsv';
Inputfile = './Dataset/Input/input_SA5_pTH-aroG_fbr-ppsA-tktA_old_190325.tsv';
%Inputfile = './Dataset/Input/input_SA5.tsv';
%Inputfile = './Dataset/Input/input_SA5_till_10h.tsv';
%Inputfile = './Dataset/Input/input_SA5_till_5h.tsv';


delimiterIn = ' ';

headerlinesIn = 1;

%A = importdata(Inputfile, delimiterIn, headerlinesIn);
A = importdata(Inputfile);

ModelFile = './Dataset/Input/iJO1366.xml';
model = readCbModel(ModelFile);

%model = changeRxnBounds(model, {'EX_o2(e)'}, [-20], 'l');

%model = changeRxnBounds(model, {'EX_o2(e)', 'EX_co2(e)'}, [-20, 0], 'l');

model = changeRxnBounds(model, {'EX_o2(e)', 'EX_co2(e)', 'EX_for(e)', 'EX_etoh(e)'}, [-20, 0, 0, 0], 'l');

%model = changeObjective(model, 'EX_skm(e)');

model_ori = model;

%n = 1201; 

[row, column] = size(A.data)

GUR = A.data(:,2);

mu = A.data(:,3);

time = A.data(:,1);

B = zeros(size(model.rxns));

size(B)


%
for i = 1: row;
	model = model_ori;
	disp(i);
	margin = 0.01; % for glc
	
	margin2 = 0.001; % for growth
	
	%Input flux values
	
	if mu(i) < 0.0
		mu(i) = 0.0;
	end;
	
	if GUR(i) >= 0.0
		GUR(i) = 0.0;
	end;
	
	if time(i) <= 2.0
		
		model = changeRxnBounds(model, {'EX_skm(e)'}, [0.0], 'u');
	end;
	
	model = changeRxnBounds(model, {'EX_glc(e)', 'Ec_biomass_iJO1366_core_53p95M'}, [GUR(i), mu(i)], 'b');
	
	solution = optimizeCbModel(model);
	
	if(solution.stat ==1)
		
		Fixed_mu = solution.x(findRxnIDs(model, 'Ec_biomass_iJO1366_core_53p95M'));
		
		model = changeObjective(model, 'EX_skm(e)');
		
		solution = optimizeCbModel(model);
		
		B = [B solution.x];
	
	else
		
		model = changeRxnBounds(model, {'EX_phe-L(e)', 'EX_tyr-L(e)', 'EX_trp-L(e)'}, [-0.5, -0.5, -0.225], 'l');
		
		solution = optimizeCbModel(model);
		
		if(solution.stat==1)
			
			Fixed_mu = solution.x(findRxnIDs(model, 'Ec_biomass_iJO1366_core_53p95M'));
			
			model = changeObjective(model, 'EX_skm(e)');
			
			solution = optimizeCbModel(model);
			
			B = [B solution.x];
			
		else
			
			for iter = 1:100;
				
				%model = changeRxnBounds(model, {'EX_glc(e)'}, [GUR(i) * (1 + margin * iter)], 'l');
				
				model = changeRxnBounds(model, {'EX_glc(e)', 'Ec_biomass_iJO1366_core_53p95M'}, [GUR(i) * (1 + margin * iter), mu(i) * (1 - margin2 * iter)], 'l');
				
				model = changeRxnBounds(model, {'EX_glc(e)'}, [GUR(i) * (1 - margin * iter)], 'u');
				
				solution = optimizeCbModel(model);
				
				if(solution.stat==1)
					
					Fixed_mu = solution.x(findRxnIDs(model, 'Ec_biomass_iJO1366_core_53p95M'));
					
					model = changeObjective(model, 'EX_skm(e)');
					
					solution = optimizeCbModel(model);
					
					B = [B solution.x];
					
					condition = ['diff = ', num2str(iter)];
					
					disp(condition);
					
					break
					
					
				end; % if end in for-iter-loop
			
			end; % for-iter-loop end
			
			if(solution.stat~=1)
				
				B = [B zeros(size(model.rxns))];
			end;
			
		end; % if sentence  end in aminoa acids fed
		
	end; % if sentence end in for-i-loop end

end; % for-i-loop end

% Output
dlmwrite('./Dataset/Output/CalcuResult_skm_SA5_pTH-aroG_fbr-ppsA-tktA_h_1e-2.tsv', B, 'delimiter', '\t');
%dlmwrite('./Dataset/Output/CalcuResult_skm_h_1e-2_till_10h.tsv', B, 'delimiter', '\t');
%dlmwrite('./Dataset/Output/CalcuResult_skm_h_1e-2_till_5h.tsv', B, 'delimiter', '\t');


