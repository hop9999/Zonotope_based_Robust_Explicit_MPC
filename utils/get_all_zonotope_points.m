function Points = get_all_zonotope_points(G)

[n, k] = size(G);

%https://www.mathworks.com/matlabcentral/answers/525295-how-do-i-create-a-matrix-with-all-binary-combinations
BinaryCubeVertices = dec2bin(0:2^k-1)' - '0';

%this is for {-1, 1} convention, without this line it would be {0, 1}
BinaryCubeVertices = 2*BinaryCubeVertices - ones(size(BinaryCubeVertices));

L = size(BinaryCubeVertices, 2);
Points = zeros(n, L);

for i = 1:L
    Points(:, i) = G*BinaryCubeVertices(:, i);
end
    
end