classdef Legendre < Recurr
    
    methods
        function obj = Legendre(order)
            k = double((0:order))';
            a = (2*k+1)./(k+1);
            b = zeros(size(k));
            c = k./(k+1);        
            normalising = reshape( sqrt(double(2*(0:order)+1)), 1, []);
            obj@Recurr(order, a, b, c, normalising);
            obj.domain = [-1,1];
            obj.constant_weight = true;
        end
        
        function x = measure_inverse_cdf(obj, z)
            x = z*2 - 1;
        end        

        %function x = sample_measure(obj, n)
        %    x = rand(1,n)*2 - 1;
        %end        

        function x = sample_measure_skip(obj, n)
            left  = (min(obj.nodes) - 1)/2;
            right = (max(obj.nodes) + 1)/2;
            x = rand(1,n)*(right-left) + left;
        end

        function w = eval_measure(obj, x)
            w = 0.5*ones(size(x));
        end        

        function w = eval_log_measure(obj, x)
            w = log(0.5)*ones(size(x));
        end   

        function w = eval_measure_deri(obj, x)
            w = zeros(size(x));
        end       
        
        function w = eval_log_measure_deri(obj, x)
            w = zeros(size(x));
        end     
    end
end