classdef CMass % By Pengyu Hong @ Brandeis University
    properties (Constant)
        cAtomNames = { 'C', 'H', 'N', 'O' };
        cAtomMasses = { [12 13.0033548378], [1.0078250321 2.0141017780], [14.0030740052, 15.0001088984], [15.9949146221 16.99913150 17.9991604]};
        cAtomProbs = { [0.9893 0.0107], [0.99985 0.00015], [0.99632, 0.00368], [0.99757 0.00038 0.00205] };
        
        Electron = 0.0005489;
        H = 1.0078250321;
        H2 = 2.0156500642;
        H2O = 18.0105646863;
        C = 12;
        N = 14.0030740052;
        O = 15.9949146221;
        CH2 = 14.0156500642;
        
        Proton = 1.007276432;
        Lithium = 7.0154553836;
        Lithium_fix = 6.8418;
        Sodium = 22.989769;
        Cesium = 132.90545;
        
        cReducingEndModification_O18 = 'O18';
        cReducingEndModification_Deuterium = 'Deuterium';
        cReducingEndModification_Reduced = 'Reduced';
        cReducingEndModification_Aminopyridine = 'Aminopyridine';
        cReducingEndModification_PRAGS = 'PRAGS';
        cReducingEndModification_M3_Bion = 'REM_M3_Bion';

        cMassCompensation_O18 = 2.00425;
        cMassCompensation_Deuterium = 17.03758;
        cMassCompensation_Aminopyridine = 78.05803471;
        cMassCompensation_PRAGS = 120.0687;
        
        cPermethylationMassLoss = CMass.CH2*2;
    end
    
    methods (Static)
        function mass = get_atom_mass( atom )
            if strcmp( atom, 'H' )
                mass = CMass.H;
            elseif strcmp( atom, 'Proton' )
                mass = CMass.Proton;
            elseif strcmp( atom, 'Na' ) || strcmp( atom, 'Sodium' )
                mass = CMass.Sodium;
            elseif strcmp( atom, 'Li' ) || strcmp( atom, 'Lithium' )
                mass = CMass.Lithium;
            elseif strcmp( atom, 'Cs' ) || strcmp( atom, 'Cesium' )
                mass = CMass.Cesium;
            elseif strcmp( atom, 'O' )
                mass = CMass.O;
            elseif strcmp( atom, 'N' )
                mass = CMass.N;
            elseif strcmp( atom, 'C' )
                mass = CMass.C;
            elseif strcmp( atom, 'Li_fix')
                mass = CMass.Lithium_fix;
            end
        end
        
        function mass = get_mass_compensation( modification, permethylated )
            if nargin < 2 || isempty( permethylated )
                permethylated = 0;
            end
            mass = 0;
            switch modification
                case CMass.cReducingEndModification_O18
                    mass = CMass.cMassCompensation_O18;
                    
                case CMass.cReducingEndModification_Deuterium
                    mass = CMass.cMassCompensation_Deuterium;
                    
                case CMass.cReducingEndModification_Reduced
                    mass = CMass.H2 + permethylated * CMass.CH2;
                    
                case CMass.cReducingEndModification_Aminopyridine
                    mass = CMass.cMassCompensation_Aminopyridine;
                
                case CMass.cReducingEndModification_PRAGS
                    mass = CMass.cMassCompensation_PRAGS;
                    
                case CMass.cReducingEndModification_M3_Bion
                    mass = -CMass.H2O;
                    if permethylated
                        mass = mass - CMass.CH2;
                    end
                    % mass = 0;
                    
                otherwise
                    if str_startswith( modification, CMass.cReducingEndModification_M3_Bion )
                        %mass = -CMass.H2O;
                        %if permethylated
                        %    mass = mass - CMass.CH2;
                        %end
                        temp = strrep( modification, [CMass.cReducingEndModification_M3_Bion, '_'], '' );
                        mass = CMass.chem_formula_2_mass( temp );
                    end
            end
        end
        
        function [monoMasses, aveMasses] = chem_formula_2_mass( formula_set )
            monoMasses = []; aveMasses = [];
            if isempty( formula_set )
                monoMasses = 0; aveMasses = 0;
                return;
            end
            
            if ischar( formula_set )
                formula_set = {formula_set};
            end
            
            for f = 1 : length( formula_set )
                formula = formula_set{f};
                
                formula = strtrim( formula );
                if isempty( formula )
                    monoMass = 0; aveMass = 0; return;
                end
                
                sign = 1;
                if formula(1) == '-'
                    sign = -1;
                    formula = formula(2:end);
                elseif formula(1) == '+'
                    formula = formula(2:end);
                end
                
                atomIdx = '';
                num = 0; monoMass = 0; aveMass = 0;
                for k = 1 : length( formula )
                    if formula(k) >= '0' && formula(k) <= '9'
                        num = num * 10 + str2num( formula(k) );
                    else
                        if ~isempty( atomIdx )
                            if num == 0, num = 1; end
                            monoMass = monoMass + num * CMass.cAtomMasses{atomIdx}(1);
                            aveMass = aveMass + num * CMass.cAtomMasses{atomIdx} * CMass.cAtomProbs{atomIdx}';
                        end
                        atomIdx = find( strcmp( CMass.cAtomNames, formula(k) ) );
                        num = 0;
                    end
                end
                
                if ~isempty( atomIdx )
                    if num == 0, num = 1; end
                    monoMass = monoMass + num * CMass.cAtomMasses{atomIdx}(1);
                    aveMass = aveMass + num * CMass.cAtomMasses{atomIdx} * CMass.cAtomProbs{atomIdx}';
                end
                monoMasses(f) = monoMass * sign;
                aveMasses(f) = aveMass * sign;
            end
        end
    end
end