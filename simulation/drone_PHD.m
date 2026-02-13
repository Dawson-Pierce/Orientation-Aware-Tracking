classdef drone_PHD < BREW.multi_target.PHD
    % The purpose of this new PHD filter is to use an additional step to
    % attempt to remove feet / extra stuff from stl file

    % NOTE: BREW must be on path to get this class to work

    methods
        function [obj,meas_new] = correct(obj, dt, meas, varargin)
                % Input is expected to be matrix of points, but added logic if
                % it's a cell of measurements
    
                % Matrix format: M x T where T is each point, and M is
                % dimension of measurement
    
                p = inputParser;
                p.KeepUnmatched = true;
                addParameter(p, 'maxDistance', 0.25);
                addParameter(p, 'referenceVector', [0 0 1]);
                addParameter(p, 'maxAngularDistance', 5);
                p.parse(varargin{:})
    
                if ~isa(meas,'cell')
                    if obj.Mix.isExtendedTarget
                        % Clustering needs to happen
                        meas_new = obj.cluster_obj.cluster(meas);
                    else
                        meas_new = mat2cell(meas, size(meas,1), ones(1, size(meas,2)));
                    end
                else
                    meas_new = meas;
                end
                % meas_new is a cell of all measurements
    
                % in the case that the target is not detected
                undetected_mix = obj.Mix.copy(); 
                undetected_mix.weights = (1 - obj.prob_detection) * undetected_mix.weights;

                % Additional point cloud manipulation for drone case 
                for k = 1:length(meas_new)
                    ptCloud_temp = pointCloud(meas_new{k}');
                    [~,inlierIndices,~] = pcfitplane(ptCloud_temp,...
                                p.Results.maxDistance,p.Results.referenceVector,p.Results.maxAngularDistance);
                    ptCloud_new = select(ptCloud_temp,inlierIndices); 
                    meas_new{k} = ptCloud_new.Location';
                end 

                % Correct the mixture for each measurement
                obj.correct_prob_density(dt, meas_new, p.Unmatched);
    
                obj.Mix.addComponents(undetected_mix);
    
        end
    end
end