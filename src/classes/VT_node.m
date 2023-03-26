classdef VT_node
    %VT_NODE Creates a node in the VT-based network representation
    
    properties
        connection;
        turningRatio;
        positions;
        offset; % regards to the first corridor mentioned in input "connection"
        cycleLength;
        redPhase;
        signal1;
        signal2=[];
    end
    
    methods
        function obj = VT_node(conn,tr, pos, os, cl, rp, simTime)
            %VT_NODE
            %   It stores which corridor are connected in this node. Also,
            %   the turning ratios are stored. Also, it makes sure the
            %   offsets are consistent.
            if (tr(2) == -1 && conn(2) ~= -1) || (tr(2) ~= -1 && conn(2) == -1)
                error('Turning ratio not well defined.');
            elseif (pos(2) == -1 && conn(2) ~= -1) || (pos(2) ~= -1 && conn(2) == -1)
                error('Node positions not well defined.');
            end
            obj.connection = conn;
            obj.turningRatio = tr;
            obj.positions = pos;
            obj.offset = os;
            obj.cycleLength = cl;
            obj.redPhase = rp;
            
            if conn(2)  ~= -1
                obj.signal1 = VT_signal(obj.cycleLength,obj.redPhase,obj.offset,obj.positions(1), simTime);
                obj.signal2 = VT_signal(obj.cycleLength,obj.redPhase,obj.offset-45,obj.positions(2), simTime);
            else
                obj.signal1 = VT_signal(obj.cycleLength,obj.redPhase,obj.offset,obj.positions(1), simTime);
            end
        end
    end
end

