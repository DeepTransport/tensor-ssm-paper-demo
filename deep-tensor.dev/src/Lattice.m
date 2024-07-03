classdef Lattice

    properties
        s
        z
    end

    methods
        function obj = Lattice(d)
            obj.s = d;
            % Getting lattice rule
            %fid = fopen('lattice-39101-1024-1048576.3600');
            %fid = fopen('lattice-38005-1024-1048576.5000');
            fid = fopen('lattice-33002-1024-1048576.9125');
            [z, ~] = fscanf(fid, '%d %d', [2, inf]);
            obj.z = z(2, 1:obj.s)';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x = generate(obj, ind, shift)
            x = mod(obj.z(:)*ind + shift(:), 1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x = new_points(obj, log2N, shift)
            if nargin == 2
                shift = zeros(obj.s,1);
            end
            if numel(shift) ~= obj.s
                error('the size of the shift vector is wrong')
            end
            n = 2^log2N;
            ind = (0:(n-1))/n;
            x = generate(obj, ind, shift);
        end

        function x = add_points(obj, log2N, shift)
            if nargin == 2
                shift = zeros(obj.s,1);
            end
            if numel(shift) ~= obj.s
                error('the size of the shift vector is wrong')
            end
            n = 2^log2N;
            ind = (1:2:(n-1))/n;
            x = generate(obj, ind, shift);
        end
    end

end
