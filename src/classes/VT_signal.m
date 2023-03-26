classdef VT_signal
    %SIGNAL Summary of this class goes here
    %   Assumption 1: Cycle always starts with a red time.
    
    properties
        cycleTime; % Duration of signal cycle
        redTime; % Duration of red time
        offset; % Duration of offset
        position; % Position
        greenTime = []; % Duration of offset
        listCycleTime=[]; % List containing the start/end times of each cycle (its the same in the case of a cycle
        listRedTime=[]; % List containing the end times of red phase
        listGreenTime=[]; % List containing the end times of green phase
    end
    
    methods
        function obj = VT_signal(cycleTime,redTime,offset,position,simulationHorizon)
            %SIGNAL Construct an instance of a signal
            if obj.redTime > obj.cycleTime
                error('Red time > cycle Time');
            elseif obj.redTime == obj.cycleTime
                error('Red time == cycle Time');
            end
            obj.cycleTime=cycleTime;
            obj.redTime=redTime; 
            obj.offset=mod(offset,cycleTime); % Makes sure only the necessary offset is considered
            obj.position=position;
            obj.greenTime = cycleTime - redTime;
            obj.listCycleTime = offset:cycleTime:simulationHorizon;
            
            % Consider of sets in the list of red and green times
            if offset > obj.greenTime
                obj.listRedTime = [redTime - (cycleTime - offset) obj.listCycleTime + redTime];
            else
                obj.listRedTime = obj.listCycleTime + redTime;
            end
            obj.listGreenTime = obj.listCycleTime;
            
        end
    end
end

