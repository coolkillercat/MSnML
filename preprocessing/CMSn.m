classdef CMSn < handle
    properties
        
        mType = '';
        mExperiment = '';
        mPrecursorMZ = -1;
        mProtonated = 0;
        mM = [];
        mPeaks = [];
        mScans = [];
        mLevel = -1; % rank of MS
        mCollisionEnergy = -1;
        mPrecursorZ = 1;
        mSumIntensity = [];
        mRawIntensityMatrix = [];
        mMixRatio = '';
        mRatio = '';
        mFeatureIntensity = [];
        %end
        mRange = [65 211];
        %properties(Access=private)
        mScanPeaks = {};
        mBatch = 1;
        mPlot = 0;
        mTemplate = [];%[89.0596214320000,6.97926414708458;103.074170717714,25.7975052681308;117.118156432000,4.57299072138346;131.070046432000,13.5156088443137;141.054736432000,3.81420775368725;143.411971432000,5.79486413592998;173.103556432000,13.7439038864885;173.273701432000,9.30152979247603;175.585906432000,6.70178224265560;175.763461432000,24.5418465053907];
        mUseComplement = false;
        mGlycan = '';
        mMetalMass = 0;
        mPrecursorMetalMass = 0;
    end

    methods
        function this = CMSn( filename, varargin)
%             if nargin == 0
%                 return
%             end
%             if nargin < 2
%                 this.load( filename );
%                 return
%             end
%             if nargin < 3
%                 this.load_with_resolution(filename, resolution);
%                 return
%             end
            p = inputParser;
            
            addRequired(p,'filename');
            addParameter(p,'resolution',100000);
            addParameter(p,'template',this.mTemplate);
            addParameter(p, 'batch', 1);
            addParameter(p, 'range', [65 211]);
            addParameter(p, 'add_complementary_peaks', false);
            addParameter(p, 'ratio', -1)
            addParameter(p, 'level', -1)
            addParameter(p, 'metal_mass', CMass.get_atom_mass('Li_fix'));
            addParameter(p, 'precursor_mz', -1);
            addParameter(p, 'precursor_metal_mass', CMass.get_atom_mass('Li'));
            parse(p,filename, varargin{:})
            this.mTemplate = p.Results.template;
            this.mBatch = p.Results.batch;
            this.mRange = p.Results.range;
            this.mRatio = p.Results.ratio;
            this.mUseComplement = p.Results.add_complementary_peaks;
            this.mLevel = p.Results.level;
            this.mMetalMass = p.Results.metal_mass;
            this.mPrecursorMZ = p.Results.precursor_mz;
            this.mPrecursorMetalMass = p.Results.precursor_metal_mass;
            if ~isempty(filename)
                this.load_with_resolution(filename, p.Results.resolution);
            end
        end
        
        % Set resolution to 10000 by default
        function load( this, filename )
            load_with_resolution( this, filename, 100000)
        end
        
        function parse_filename(this, filename)
            idx = strfind(filename,'_rerun');
            if isempty(idx)
                idx = strfind(filename, '.mzXML');
            end
            if contains(filename, 'CID')
                this.mExperiment = 'CID';
                idxs = strfind(filename, 'CID');
                this.mType = filename(idxs-1);
                this.mMixRatio = filename((idxs+6):(idx-1));
            end
            if contains(filename, 'HCD')
                this.mExperiment = 'HCD';
                idxs = strfind(filename, 'HCD');
                this.mType = filename(idxs-1);
                this.mMixRatio = filename((idxs+6):(idx-1));
            end
        end
        
        function temp = init(this, filename)
            this.mM = [];
            this.mScans = [];
            this.mScanPeaks = {};
            this.mPeaks = {};
            this.mProtonated = 0;
            temp = mzxmlread( filename );
            this.parse_filename(filename);
            if this.mPrecursorMZ == -1
                this.mPrecursorMZ = temp.scan(1).precursorMz.value;
                this.mPrecursorZ = temp.scan(1).precursorMz.precursorCharge;
            end
            if this.mLevel == -1
                this.mLevel = temp.scan(1).msLevel;
            end
            this.mCollisionEnergy = temp.scan(1).collisionEnergy;
        end
        
        function [X, Intensities] = resample_peak(this, temp, resolution)
            [Peaklist, ~] = mzxml2peaks(temp,'Levels',this.mLevel);

            if this.mUseComplement == true
                list_l = size(Peaklist, 1);
                for i = 1:list_l
                    curr = Peaklist{i};
                    l1 = this.mPrecursorMZ - curr(end:-1:1, 1) + this.mMetalMass;
                    l2 = curr(end:-1:1, 2);
                    Peaklist{end+1} = [l1 l2];
                end
            end
            [X, Intensities] = msppresample(Peaklist, resolution, 'Range', this.mRange,'FWHH',0.25);
            Intensities(1,:)=zeros(1, size(Intensities,2));
            Intensities(end,:)=zeros(1, size(Intensities,2));
            %[X, Intensities] = msresample(X, Intensities, resolution,'Uniform', 'true', 'Range', this.mRange);
            Intensities = max(0, Intensities);
            this.mRawIntensityMatrix = Intensities;
        end
        
        function Intensities_batch = batch_intensities(this, Intensities, batch)
            Intensities_batch = [];
            for i = 1:floor(size(Intensities,2)/batch)
                subIntensities = Intensities(:, ((i-1)*batch+1):(i*batch));
                Intensities_batch(:,end+1) = sum(subIntensities,2);
            end
        end
        
        function add_complementary_peaks(this, X, Intensities)
            
        end
        
        function load_with_resolution( this, filename, resolution)
            temp = this.init(filename);
            [X, Intensities] = this.resample_peak(temp, resolution);
            if this.mUseComplement == true
                this.add_complementary_peaks(X, Intensities);
            end
            this.mM = X;
            this.protonate(this.mPrecursorMetalMass, this.mMetalMass);
            X = this.mM;
            %Intensities(Intensities < 0.1) = 0;
            Intensities = batch_intensities(this, Intensities, this.mBatch);
            Intensities = Intensities(:, sum(Intensities,1) > 0);
            if (~isempty(this.mTemplate))
            YA = msalign(X, Intensities, this.mTemplate(:,1), 'ShowPlot', 'false', 'MaxShift', [-0.5 0.5], 'WidthOfPulses', 0.05, 'WEIGHTS', this.mTemplate(:,2));
            Intensities = YA;
            end  
            Intensities(isnan(Intensities)) = 0; 
            %Intensities = msalign(X, Intensities, reference, 'ShowPlot','true', 'MaxShift', [-1 1], 'WidthOfPulses', 0.1);
            if isempty( temp ) || isempty( temp.scan )
                return;
            end
            
            this.mM = X;
            this.mScans = Intensities;
            temp = sum( this.mScans, 2 );
            this.mSumIntensity = sum(this.mScans, 2);
        end
        
        
        
        function detect_peaks( this )
            temp = this.mSumIntensity;%sum( this.mScans, 2 );
            th = size(this.mScans, 2)/5;
            if this.mLevel > 2
                this.mPeaks = mspeaks( this.mM, temp, 'HeightFilter', 5 );
            else
                temp = mspeaks( this.mM, temp, 'HeightFilter', 1 );
                N = size(temp, 1);
                flag = false(N, 1);
                iso_flag = false(N, 1);
                for m = 1 : N-1
                    if temp(m, 1) < th
                        continue;
                    end
                    iso = temp(m, 1) + 1;
                    [mv, idx] = min(abs(temp(:,1) - iso));
                    if mv < 0.0667 && temp(m,2) > temp(idx,2) * 1.25
                        iso2 = temp(m, 1) + 2;
                        iso_flag(idx) = true;
                        [mv2, idx2] = min(abs(temp(:,1) - iso2));
                        if mv2 < 0.134 % && abs(temp(idx2,1) - temp(idx, 1) - 1) < 0.07
                            flag(m) = true;
                            if temp(idx,2) > temp(idx2, 2) * 1.25
                                iso_flag(idx2) = true;
                            end
                        end
                    end
                end
                this.mPeaks = temp(flag & ~iso_flag, :);
            end
        end
        
        function normalizePeaks(this)
            if isempty(this.mPeaks)
                return
            end
%             maxIntensity = max(this.mPeaks(:,2));
%             for i = 1:length(this.mPeaks)
%                 this.mPeaks(i,2) = this.mPeaks(i,2) / maxIntensity;
%             end
            this.mPeaks(:,2) = msnorm(this.mPeaks(:,1),this.mPeaks(:,2),'Quantile',[0.75 1],'Consensus',0.9);
%             maxIntensity = max(this.mPeaks(:,2));
%             for i = 1:length(this.mPeaks)
%                 this.mPeaks(i,2) = this.mPeaks(i,2) / maxIntensity;
%             end
%             P = prctile(this.mPeaks(:,2), 0.95);
%             term = find(this.mPeaks(:,2) >= P);
%             sumOfI = sum(this.mPeaks(term, 2));
%             for i = 1:length(this.mPeaks)
%                 this.mPeaks(i,2) = this.mPeaks(i,2) / sumOfI;
%             end
        end

        function aPeak = get_peak( this, mass, massAccuracy )
            [mv, idx] = min( abs( this.mPeaks(:,1) - mass ) );
            if mv < massAccuracy
                aPeak = this.mPeaks(idx,:);
            else
                aPeak = [];
            end
        end

        function align_scanpeaks( this )
            numScan = size(this.mScans, 2);
            if isempty( this.mScanPeaks )
                this.detect_scanpeaks();
            end
            template = this.mScanPeaks{1};
            count = ones(size(template,1), 1);
            for k = 2 : numScan
                [template, ~, count] = CMSn.align_two_peaklists( template, this.mScanPeaks{k}, 0.2, 10, count );
            end
            this.mPeaks = template;
        end

        function detect_scanpeaks( this )
            numScan = size(this.mScans, 2);
            mask = sum( this.mScans > 0, 2 ) > 0.2 * numScan;
            this.mScanPeaks = cell(numScan, 1);
            for k = 1 : numScan
                temp = mspeaks( this.mM, this.mScans(:, k) .* mask, 'HeightFilter', 1 );
                N = size(temp, 1);
                flag = false(N, 1);
                for m = 1 : N-1
                    iso = temp(m, 1) + 1;
                    [mv, idx] = min(abs(temp(:,1) - iso));
                end
                this.mScanPeaks{k} = temp(flag>0, :);
            end
        end

        function protonate_precursor( this, ionmass, precursorZ) % for example, ionmass = CMass.get_atom_mass("Li")
            if nargin < 3, precursorZ = this.mPrecursorZ; end
            if this.mProtonated == 0
                this.mPrecursorMZ = (this.mPrecursorMZ - ionmass) * precursorZ + CMass.Proton;
                this.mPrecursorZ = 1;
            end
        end

        function protonate_signal( this, ionmass ) % for example, ionmass = CMass.Lithium - 0.1736
            if this.mProtonated == 1
                return
            end
            this.mM = this.mM - ionmass + CMass.Proton;
            if ~isempty( this.mPeaks )
                this.mPeaks(:,1) = this.mPeaks(:,1) - ionmass + CMass.Proton;
            end
        end

        function protonate(this, precursorIonmass, peakIonmass, precursorZ)
            if nargin < 4
                this.protonate_precursor(precursorIonmass);
            else
                this.protonate_precursor(precursorIonmass, precursorZ);
            end
            this.protonate_signal(peakIonmass);
            this.mProtonated = 1;
        end

        function scanvalues = get_intensities( this, mz )
            [d, idx] = min( abs(this.mM - mz) );
            if (d < 0.1)
                scanvalues = this.mScans(idx, :);
            else
                scanvalues = [];
            end
        end

        function visualize_scans( this, scanIdxes )
            figure;
            if nargin < 2 || isempty( scanIdxes )
                stem( this.mM, this.mSumIntensity   , 'marker', 'none' );
            else
                stem( this.mM, this.mScans(:, scanIdxes), 'marker', 'none' );
            end
        end
        
        function visualize_scanpeaks( this, scanIdx )
            figure;
            stem( this.mScanPeaks{scanIdx}(:,1), this.mScanPeaks{scanIdx}(:,2), 'marker', 'none' );
        end

        function visualize_peaks( this )
            figure;
            stem( this.mPeaks(:,1), this.mPeaks(:,2), 'marker', 'none' );
        end
        
        function heatmap( this )
            msheatmap(this.mM, 1:size(this.mScans, 2), log(this.mScans+1));
            ylabel( 'Scan ID' );
        end
        
        function f = check_features(this, lo, up)
            f = [];
            if (isempty(this.mPeaks))
                f = zeros(length(lo),1);
                return
            end
            for i = 1:length(lo)
                m = this.mPeaks(:,1);
                m = m(m >= lo(i));
                m = m(m <= up(i));
                if ~isempty(m)
                    f(end+1) = 1;
                else
                    f(end+1) = 0;
                end
            end
        end
        
        function f = check_features_by_threshold(this, features, threshold)
            f = zeros(1,length(features));
            if (isempty(this.mPeaks) || isempty(features))
                return
            end
            %normalized = normalize(log(this.mPeaks(:, 2) + 1),'range',[0 1])s;
            for i = 1:size(this.mPeaks, 1)
                m = this.mPeaks(i,1);
                [minv, idx] = min(abs(m - features));
                if minv < threshold
                    f(idx) = f(idx) + this.mPeaks(i,2);
                end
%                     idxes = find(abs(m - features) < threshold/2);
%                     for idx = idxes
%                         f(idx) = f(idx) + normpdf(abs(m - features(idx))/threshold) * this.mPeaks(i, 2);
%                     end

            end
            this.mFeatureIntensity = f;
        end
        
       function f = check_features_by_threshold_min(this, features, threshold)
            f = zeros(1,length(features));
            if (isempty(this.mPeaks) || isempty(features))
                return
            end
            %normalized = normalize(log(this.mPeaks(:, 2) + 1),'range',[0 1]);
            for i = 1:size(this.mPeaks, 1)
                m = this.mPeaks(i,1);
                 [minv, idx] = min(abs(m - features));
                 if minv < threshold
%                     idxes = find(abs(m - features) < 3 * threshold);
%                     for idx = idxes
%                         f(idx) = f(idx) + normpdf(abs(m - features(idx))/threshold) * this.mPeaks(i, 2);
%                     end
                    f(idx) = f(idx) + this.mPeaks(i, 2);
                 end
            end
            this.mFeatureIntensity = f;
       end
        
       function apply_medfilter(this, n)
           if (isempty(this.mSumIntensity)) return; end
           this.mSumIntensity = medfilt1(this.mSumIntensity, n);
           return
       end
        
    end

    methods (Static)
        function [alignedSpec, overlaps, count] = align_two_peaklists( spec1, spec2, massAccuracy, intensityStd, count )
        % This function is not very good yet.
            if nargin < 3 || isempty(massAccuracy)
                massAccuracy = 0.01;
            end
            if nargin < 4 || isempty(intensityStd)
                intensityStd = 10000;
            end

            n1 = size( spec1, 1 );
            n2 = size( spec2, 1 );
            M = zeros( n1+1, n2+1 );

            if nargin < 5 || isempty(count)
                count = ones(n1, 1);
            end

            M(n1+1, :) = 0.001;
            M(:, n2+1) = 0.001;
            for k = 1 : n1
                for m = 1 : n2
                    dMass = spec1(k, 1) - spec2(m, 1);
                    if dMass > massAccuracy || dMass < -massAccuracy
                        continue;
                    end
                    dIntensity = abs(spec1(k, 2) - spec2(m, 2));
                    M(k, m) = exp(-dMass/massAccuracy) * exp( -dIntensity / intensityStd );
                end
            end

            mask = M > 0;
            for beta = 5 : 0.5 : 10
                Q = exp(beta * M) .* mask;
                QQ = Q;
                for k = 1 : n1
                    Q(k,:) = Q(k,:) / sum( Q(k,:));
                end
                for m = 1 : n2
                    QQ(:,m) = QQ(:,m) / sum( QQ(:,m));
                end
                newQ = (Q + QQ) / 2;
                newQ(n1+1, :) = QQ(n1+1, :);
                newQ(:, n2+1) = Q(:, n2+1);
                d = max(max(abs(newQ-M)));
                M = newQ;
                if d < 1e-6
                    break
                end
            end

            alignedSpec = spec1;
            matching = zeros(n1, 1);
            matched = zeros(n2, 1);
            for k = 1 : n1
                [v, idx] = max(M(k,:));
                if idx <= n2 && v > 0.5
                    matching(k) = idx;
                    mass = (spec1(k,1) * count(k) + spec2(idx,1)) / (count(k) + 1);
                    alignedSpec(k,1) = mass;
                    alignedSpec(k,2) = (alignedSpec(k,2) * count(k) + spec2(idx,2)) / (count(k) + 1);
                    matched(idx) = 1;
                    count(k) = count(k) + 1;
                end
            end
            alignedSpec = [alignedSpec; spec2(matched == 0, :)];
            numNew = sum(matched == 0);
            overlaps = [matching > 0; zeros(numNew, 1)];
            count = [count; ones(numNew,1)];
            [alignedSpec, idxes] = sortrows(alignedSpec);
            overlaps = overlaps(idxes, :);
            count = count(idxes, :);
        end

        function mergedPeaks = merge_peaks( spec, massTolerance )
            if nargin < 2 || isempty( massTolerance )
                massTolerance = 0.01;
            end
            while true
                d = spec(2:end, 1) - spec(1:end-1, 1);
                [mv, idx] = min(d);
                if mv < massTolerance
                    newIntensity = spec(idx,2) + spec(idx+1,2);
                    newMass = (spec(idx,1)*spec(idx,2) + spec(idx+1,1)*spec(idx+1,2)) / newIntensity;
                    spec(idx, :) = [newMass, newIntensity];
                    spec = [spec(1:idx, :); spec(idx+2:end, :)];
                else
                    break;
                end
            end
            mergedPeaks = spec;
        end
        
        function new = filter(old, type, method, power, mix)
            new = {};
            for i = 1:length(old)
                MS = old{i};
                if(~isempty(type) && MS.mType ~= type) continue; end
                if(~isempty(method) && strcmp(MS.mExperiment, method) == 0) continue; end
                if(~isempty(power) && ~any(MS.mCollisionEnergy == power)) continue; end
                if(~isempty(mix) && strcmp(MS.mMixRatio, mix) == 0) continue; end
                new{end+1} = MS;
            end
        end
        
        function quickplot(MS5s)
            hold on
            for i = 1:length(MS5s)
                plot(MS5s{i}.mM, MS5s{i}.mSumIntensity);
            end
            hold off
        end
        
        function quickplot_raw(MS5s)
            hold on
            for i = 1:length(MS5s)
                plot(MS5s{i}.mM, sum(MS5s{i}.mRawIntensityMatrix,2));
            end
            hold off
        end
       
        function quickplot_relative(MS5s)
            sum = zeros(length(MS5s{1}.mM),1);
            for i = 1:length(MS5s)
                if max(MS5s{i}.mSumIntensity) == 0 
                    continue
                end
                sum = sum + MS5s{i}.mSumIntensity / max(MS5s{i}.mSumIntensity);
            end
            plot(MS5s{1}.mM, sum)
        end
        
        function [SL3, SL9to1, SL1to1, SL1to9, SL6] = group_by_ratio(MS5s)
            SL3 = CMSn.filter(MS5s,[],[],[],'3SL');
            SL9to1 = CMSn.filter(MS5s,[],[],[],'9to1');
            SL1to1 = CMSn.filter(MS5s,[],[],[],'1to1');
            SL1to9 = CMSn.filter(MS5s,[],[],[],'1to9');
            SL6 = CMSn.filter(MS5s,[],[],[],'6SL');
        end
        
        function new = filter_idx(old, type, method, power, mix)
            new = [];
            for i = 1:length(old)
                MS = old{i};
                if(~isempty(type) && MS.mType ~= type) continue; end
                if(~isempty(method) && strcmp(MS.mExperiment, method) == 0) continue; end
                if(~isempty(power) && ~any(MS.mCollisionEnergy == power)) continue; end
                if(~isempty(mix) && strcmp(MS.mMixRatio, mix) == 0) continue; end
                new(end+1) = i;
            end
        end
        
                
        function new = get_glycan(old, gtype)
            new = {};
            for i = 1:length(old)
                MS = old{i};
                if(~isempty(gtype) && strcmp(MS.mGlycan,gtype) == 0) continue; end
                new{end+1} = MS;
            end
        end
    end
end
