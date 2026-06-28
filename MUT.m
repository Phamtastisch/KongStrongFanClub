classdef MUT

    properties
        sobj_struct (1,1) sparameters
        S11_param {mustBeNumeric}
        S21_param {mustBeNumeric}
    end

    methods
        function obj = MUT(sobj)
            obj.sobj_struct = sobj;
            obj.S11_param = squeeze(sobj.Parameters(1,1,:));
            obj.S21_param = squeeze(sobj.Parameters(1,2,:));
        end
    end

end    