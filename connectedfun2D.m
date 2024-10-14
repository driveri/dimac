% 2D Connectivity analysis
% Ideally the inputted map should be a binary mask
%
% b = connectedfun2D(a)
%

function b = connectedfun2D(a)
s = size(a);
b = zeros(size(a));

n = 1;
%seedpos = zeros(2,240*240*62);
SE = ones(3,3);  % 2D structure element



    for y = 1:s(2)

        for x = 1:s(1)

		if a(x,y) ~= 0 && b(x,y) == 0
			seed = zeros(size(a));
			seed(x,y) =1;  % seed voxel to start the iteration

			connected_iteration = seed;
			convergence_test = zeros(size(a));

			while sum(connected_iteration(:)) ~= sum(convergence_test(:))

				convergence_test = connected_iteration;

				connected_iteration = double((conv2(connected_iteration,SE,'same').*a)>0);

			end

			b = b + connected_iteration*n;
%			seedpos(:,n) = [x,y];
			n = n+1;
		end

        end
    if false % set true for debugging
        disp(['Column number ',num2str(y),' out of ',num2str(s(2)),' complete.'])
    end
    end
