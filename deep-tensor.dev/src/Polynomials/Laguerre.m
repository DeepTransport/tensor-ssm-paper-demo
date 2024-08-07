classdef Laguerre < Recurr
    
    methods
        function obj = Laguerre(order)
            k = double((0:order))';
            a = -1./(k+1);
            b = (2*k+1)./(k+1);
            c = k./(k+1);
            normalising = ones(1, order+1);
            obj@Recurr(order, a, b ,c, normalising);
            obj.domain = [0,inf];
            obj.constant_weight = false;
            %
        end
        
        function x = measure_inverse_cdf(obj, z)
            tmp = 1-z;
            tmp(tmp<1E-10) = 1E-10;
            tmp(tmp>1-1E-10) = 1-1E-10;
            x = -log(tmp);
        end        

        %function x = sample_measure(obj, n)
        %    x = exprnd(1,1,n);
        %end      
  
        function w = eval_measure(obj, x)
            w = exp(-x);
        end  

        function x = sample_measure_skip(obj, n)
            x = sample_measure(obj, n);
        end
        
        function w = eval_log_measure(obj, x)
            w = -x;
        end   

        function w = eval_measure_deri(obj, x)
            w = -exp(-x);
        end       
        
        function w = eval_log_measure_deri(obj, x)
            w = -ones(size(x));
        end     
    end
    
end